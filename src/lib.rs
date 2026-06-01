use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

/// A single-cell count matrix in 10x MatrixMarket layout: rows are genes,
/// columns are cells, stored as coordinate triplets. Counts are held as f64
/// because scanpy normalizes after promoting the integer matrix to float.
pub struct CountMatrix {
    pub n_genes: usize,
    pub n_cells: usize,
    /// One entry per stored nonzero, in the file's row-major (gene-major) order.
    pub entries: Vec<Entry>,
}

#[derive(Clone, Copy)]
pub struct Entry {
    pub gene: u32,
    pub cell: u32,
    pub value: f64,
}

pub struct NormalizeParams {
    /// `None` reproduces scanpy's default: the median of every cell's total
    /// count, computed over all cells including those with zero counts.
    pub target_sum: Option<f64>,
    pub log1p: bool,
}

/// Resolve the 10x triple inside `dir`, accepting both gzipped (v3) and plain
/// (v2) layouts. Only the matrix is required for normalization.
pub fn open_mtx(dir: &Path) -> Result<Box<dyn Read>> {
    for name in ["matrix.mtx.gz", "matrix.mtx"] {
        let path = dir.join(name);
        if path.exists() {
            return open_maybe_gz(&path);
        }
    }
    Err(RsomicsError::InvalidInput(format!(
        "no matrix.mtx or matrix.mtx.gz in {}",
        dir.display()
    )))
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Parse a MatrixMarket coordinate file (real, integer, or pattern; general).
/// 10x stores genes on rows, cells on columns.
pub fn parse_mtx(reader: impl Read) -> Result<CountMatrix> {
    let mut reader = BufReader::new(reader);
    let mut line = String::new();

    reader.read_line(&mut line).map_err(RsomicsError::Io)?;
    let banner = line.trim();
    if !banner.starts_with("%%MatrixMarket") {
        return Err(RsomicsError::InvalidInput(
            "missing %%MatrixMarket banner".into(),
        ));
    }
    let pattern = banner.contains("pattern");

    let (n_genes, n_cells, nnz) = loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(RsomicsError::Io)?;
        if n == 0 {
            return Err(RsomicsError::InvalidInput("truncated MTX header".into()));
        }
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }
        let mut it = t.split_whitespace();
        let rows = parse_usize(it.next())?;
        let cols = parse_usize(it.next())?;
        let nnz = parse_usize(it.next())?;
        break (rows, cols, nnz);
    };

    let mut entries = Vec::with_capacity(nnz);
    for raw in reader.lines() {
        let raw = raw.map_err(RsomicsError::Io)?;
        let t = raw.trim();
        if t.is_empty() {
            continue;
        }
        let mut it = t.split_whitespace();
        let gene = parse_usize(it.next())?;
        let cell = parse_usize(it.next())?;
        let value = if pattern {
            1.0
        } else {
            it.next()
                .ok_or_else(|| RsomicsError::InvalidInput("MTX entry missing value".into()))?
                .parse::<f64>()?
        };
        if gene == 0 || gene > n_genes || cell == 0 || cell > n_cells {
            return Err(RsomicsError::InvalidInput(format!(
                "MTX index out of bounds: ({gene}, {cell})"
            )));
        }
        entries.push(Entry {
            gene: (gene - 1) as u32,
            cell: (cell - 1) as u32,
            value,
        });
    }
    if entries.len() != nnz {
        return Err(RsomicsError::InvalidInput(format!(
            "MTX declared {nnz} entries, found {}",
            entries.len()
        )));
    }

    Ok(CountMatrix {
        n_genes,
        n_cells,
        entries,
    })
}

/// scanpy's `np.median` over the per-cell totals (linear interpolation: for an
/// even count the two central order statistics are averaged). Zero-count cells
/// are part of the population.
fn median(totals: &[f64]) -> f64 {
    let mut sorted: Vec<f64> = totals.to_vec();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    }
}

/// Per-cell total counts, indexed by cell.
fn cell_totals(m: &CountMatrix) -> Vec<f64> {
    let mut totals = vec![0.0_f64; m.n_cells];
    for e in &m.entries {
        totals[e.cell as usize] += e.value;
    }
    totals
}

/// Normalize each cell to `target_sum` (or the median of totals) then optionally
/// apply `ln(1+x)`. Mirrors scanpy: a zero-count cell's scaling factor collapses
/// to 1, leaving its (empty) row untouched, and the sparsity pattern is exact
/// because `log1p(0) = 0`.
pub fn normalize(m: &mut CountMatrix, params: &NormalizeParams) {
    let totals = cell_totals(m);
    let target = params.target_sum.unwrap_or_else(|| median(&totals));

    let scale: Vec<f64> = totals
        .iter()
        .map(|&t| {
            let s = t / target;
            if s == 0.0 { 1.0 } else { s }
        })
        .collect();

    let log1p = params.log1p;
    m.entries.par_iter_mut().for_each(|e| {
        let v = e.value / scale[e.cell as usize];
        e.value = if log1p { v.ln_1p() } else { v };
    });
}

/// Write the matrix back in genes×cells MatrixMarket real coordinate layout,
/// preserving the input entry order. A big buffer plus ryu float formatting
/// keeps the matrix-sized write I/O-bound rather than format-bound.
pub fn write_mtx(m: &CountMatrix, out: impl Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    w.write_all(b"%%MatrixMarket matrix coordinate real general\n")
        .map_err(RsomicsError::Io)?;
    let mut header = itoa_line(m.n_genes, m.n_cells, m.entries.len());
    header.push('\n');
    w.write_all(header.as_bytes()).map_err(RsomicsError::Io)?;

    let mut fmt = ryu::Buffer::new();
    let mut buf: Vec<u8> = Vec::with_capacity(64);
    for e in &m.entries {
        buf.clear();
        write_uint(&mut buf, e.gene as u64 + 1);
        buf.push(b' ');
        write_uint(&mut buf, e.cell as u64 + 1);
        buf.push(b' ');
        buf.extend_from_slice(fmt.format(e.value).as_bytes());
        buf.push(b'\n');
        w.write_all(&buf).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

fn itoa_line(a: usize, b: usize, c: usize) -> String {
    format!("{a} {b} {c}")
}

fn write_uint(buf: &mut Vec<u8>, mut n: u64) {
    if n == 0 {
        buf.push(b'0');
        return;
    }
    let start = buf.len();
    while n > 0 {
        buf.push(b'0' + (n % 10) as u8);
        n /= 10;
    }
    buf[start..].reverse();
}

fn parse_usize(tok: Option<&str>) -> Result<usize> {
    tok.ok_or_else(|| RsomicsError::InvalidInput("MTX header missing a dimension".into()))?
        .parse::<usize>()
        .map_err(Into::into)
}

/// End-to-end: read the 10x matrix from `dir`, normalize, write to `out`.
pub fn run(dir: &Path, params: &NormalizeParams, out: impl Write) -> Result<(usize, usize)> {
    let mut m = parse_mtx(open_mtx(dir)?)?;
    let shape = (m.n_genes, m.n_cells);
    normalize(&mut m, params);
    write_mtx(&m, out)?;
    Ok(shape)
}

/// `--target-sum` accepts a positive float or the literal `median`.
pub fn parse_target_sum(s: &str) -> Result<Option<f64>> {
    if s.eq_ignore_ascii_case("median") {
        return Ok(None);
    }
    let v = s
        .parse::<f64>()
        .map_err(|_| RsomicsError::InvalidInput(format!("invalid --target-sum '{s}'")))?;
    if v <= 0.0 || !v.is_finite() {
        return Err(RsomicsError::InvalidInput(
            "--target-sum must be a positive finite number or 'median'".into(),
        ));
    }
    Ok(Some(v))
}

/// Output destination — stdout for `-`, otherwise a file.
pub fn open_output(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" {
        Ok(Box::new(std::io::stdout().lock()))
    } else {
        Ok(Box::new(
            File::create(PathBuf::from(path)).map_err(RsomicsError::Io)?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny() -> CountMatrix {
        let mut entries = Vec::new();
        let push = |v: &mut Vec<Entry>, g: u32, c: u32, val: f64| {
            v.push(Entry {
                gene: g,
                cell: c,
                value: val,
            })
        };
        push(&mut entries, 0, 0, 3.0);
        push(&mut entries, 2, 0, 1.0);
        push(&mut entries, 4, 0, 2.0);
        push(&mut entries, 1, 1, 5.0);
        push(&mut entries, 0, 2, 1.0);
        push(&mut entries, 1, 2, 1.0);
        push(&mut entries, 2, 2, 1.0);
        push(&mut entries, 3, 2, 1.0);
        CountMatrix {
            n_genes: 5,
            n_cells: 4,
            entries,
        }
    }

    #[test]
    fn median_includes_zero_cells() {
        assert_eq!(median(&[6.0, 5.0, 4.0, 0.0]), 4.5);
    }

    #[test]
    fn matches_scanpy_tiny() {
        let mut m = tiny();
        normalize(
            &mut m,
            &NormalizeParams {
                target_sum: None,
                log1p: true,
            },
        );
        let want = [
            (0u32, 0u32, 1.178655_f64),
            (2, 0, 0.5596158),
            (4, 0, 0.91629076),
            (1, 1, 1.7047482),
            (0, 2, 0.7537718),
            (1, 2, 0.7537718),
            (2, 2, 0.7537718),
            (3, 2, 0.7537718),
        ];
        for (e, (_, _, exp)) in m.entries.iter().zip(want.iter()) {
            assert!((e.value - exp).abs() < 1e-5, "{} vs {}", e.value, exp);
        }
    }

    #[test]
    fn target_sum_parsing() {
        assert_eq!(parse_target_sum("median").unwrap(), None);
        assert_eq!(parse_target_sum("1e4").unwrap(), Some(10000.0));
        assert!(parse_target_sum("-1").is_err());
        assert!(parse_target_sum("abc").is_err());
    }

    #[test]
    fn roundtrip_mtx() {
        let m = tiny();
        let mut buf = Vec::new();
        write_mtx(&m, &mut buf).unwrap();
        let parsed = parse_mtx(&buf[..]).unwrap();
        assert_eq!(parsed.n_genes, 5);
        assert_eq!(parsed.n_cells, 4);
        assert_eq!(parsed.entries.len(), 8);
    }
}

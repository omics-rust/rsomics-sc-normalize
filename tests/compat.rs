use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const EPSILON: f64 = 1e-5;

fn scanpy_python() -> Option<PathBuf> {
    let home = std::env::var("HOME").ok()?;
    let shared = PathBuf::from(&home).join("oracle-venvs/scanpy/bin/python");
    if shared.exists() {
        return Some(shared);
    }
    // fall back to a python that can import scanpy
    for cand in ["python3", "python"] {
        let ok = Command::new(cand)
            .args(["-c", "import scanpy"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if ok {
            return Some(PathBuf::from(cand));
        }
    }
    None
}

fn parse_mtx(text: &str) -> HashMap<(u64, u64), f64> {
    let mut map = HashMap::new();
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let banner = lines.next().unwrap();
    assert!(banner.starts_with("%%MatrixMarket"), "bad banner: {banner}");
    let mut dim_seen = false;
    for line in lines {
        let t = line.trim();
        if t.starts_with('%') {
            continue;
        }
        let mut it = t.split_whitespace();
        let a: u64 = it.next().unwrap().parse().unwrap();
        let b: u64 = it.next().unwrap().parse().unwrap();
        match it.next() {
            Some(v) if dim_seen => {
                map.insert((a, b), v.parse::<f64>().unwrap());
            }
            _ => {
                // dimension line (rows cols nnz)
                dim_seen = true;
            }
        }
    }
    map
}

#[test]
fn matches_scanpy_value_level() {
    let Some(py) = scanpy_python() else {
        eprintln!("SKIP: scanpy venv not found (~/oracle-venvs/scanpy/bin/python); compat skipped");
        return;
    };

    let manifest = env!("CARGO_MANIFEST_DIR");
    let mtx_dir = Path::new(manifest).join("tests/golden/tenx");
    let oracle_py = Path::new(manifest).join("tests/scanpy_normalize_oracle.py");
    assert!(mtx_dir.exists(), "missing golden 10x dir {mtx_dir:?}");

    let scratch = std::env::var("RSOMICS_SCRATCH")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::env::temp_dir());
    let oracle_out = scratch.join("sc_norm_oracle.mtx");

    let status = Command::new(&py)
        .arg(&oracle_py)
        .arg(&mtx_dir)
        .arg(&oracle_out)
        .arg("median")
        .status()
        .expect("run scanpy oracle");
    assert!(status.success(), "oracle failed");

    let ours = Command::new(env!("CARGO_BIN_EXE_rsomics-sc-normalize"))
        .arg(&mtx_dir)
        .arg("-o")
        .arg("-")
        .arg("-q")
        .output()
        .expect("run rsomics-sc-normalize");
    assert!(
        ours.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let oracle_map = parse_mtx(&std::fs::read_to_string(&oracle_out).unwrap());
    let ours_map = parse_mtx(&String::from_utf8(ours.stdout).unwrap());

    // Sparsity pattern must be exact.
    assert_eq!(
        oracle_map.len(),
        ours_map.len(),
        "nnz differs: oracle {} vs ours {}",
        oracle_map.len(),
        ours_map.len()
    );
    let mut max_dev = 0.0_f64;
    for (k, &ov) in &oracle_map {
        let our = ours_map
            .get(k)
            .unwrap_or_else(|| panic!("entry {k:?} present in oracle, missing in ours"));
        let dev = (ov - our).abs();
        max_dev = max_dev.max(dev);
        assert!(dev < EPSILON, "value mismatch at {k:?}: {ov} vs {our}");
    }
    for k in ours_map.keys() {
        assert!(
            oracle_map.contains_key(k),
            "entry {k:?} present in ours, missing in oracle"
        );
    }
    eprintln!(
        "compat OK: {} nonzeros, max deviation {max_dev:e}",
        ours_map.len()
    );
}

// Always-run guard against a committed scanpy golden (normalize_total median +
// log1p, captured from scanpy 1.11.5). Keeps CI honest where the venv is absent.
#[test]
fn matches_committed_golden() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let mtx_dir = Path::new(manifest).join("tests/golden/tenx");
    let golden = Path::new(manifest).join("tests/golden/scanpy_normalize.mtx");

    let ours = Command::new(env!("CARGO_BIN_EXE_rsomics-sc-normalize"))
        .arg(&mtx_dir)
        .arg("-o")
        .arg("-")
        .arg("-q")
        .output()
        .expect("run rsomics-sc-normalize");
    assert!(
        ours.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let golden_map = parse_mtx(&std::fs::read_to_string(&golden).unwrap());
    let ours_map = parse_mtx(&String::from_utf8(ours.stdout).unwrap());

    assert_eq!(
        golden_map.len(),
        ours_map.len(),
        "nnz differs: golden {} vs ours {}",
        golden_map.len(),
        ours_map.len()
    );
    let mut max_dev = 0.0_f64;
    for (k, &gv) in &golden_map {
        let our = ours_map
            .get(k)
            .unwrap_or_else(|| panic!("entry {k:?} in golden, missing in ours"));
        let dev = (gv - our).abs();
        max_dev = max_dev.max(dev);
        assert!(dev < EPSILON, "value mismatch at {k:?}: {gv} vs {our}");
    }
    eprintln!(
        "golden OK: {} nonzeros, max deviation {max_dev:e}",
        ours_map.len()
    );
}

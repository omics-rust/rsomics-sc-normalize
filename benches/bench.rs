use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_sc_normalize::{CountMatrix, Entry, NormalizeParams, normalize};

fn synthetic(n_genes: usize, n_cells: usize, density: f64) -> CountMatrix {
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut entries = Vec::new();
    for cell in 0..n_cells {
        for gene in 0..n_genes {
            if (next() % 10_000) as f64 / 10_000.0 < density {
                let v = (next() % 50 + 1) as f64;
                entries.push(Entry {
                    gene: gene as u32,
                    cell: cell as u32,
                    value: v,
                });
            }
        }
    }
    CountMatrix {
        n_genes,
        n_cells,
        entries,
    }
}

fn bench_normalize(c: &mut Criterion) {
    let m = synthetic(2000, 5000, 0.05);
    c.bench_function("normalize_log1p_2000x5000_5pct", |b| {
        b.iter_batched(
            || CountMatrix {
                n_genes: m.n_genes,
                n_cells: m.n_cells,
                entries: m.entries.clone(),
            },
            |mut m| {
                normalize(
                    &mut m,
                    &NormalizeParams {
                        target_sum: None,
                        log1p: true,
                    },
                );
                m
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, bench_normalize);
criterion_main!(benches);

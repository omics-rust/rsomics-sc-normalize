# rsomics-sc-normalize

Library-size normalization plus `log1p` of a single-cell count matrix,
numerically matching scanpy's `sc.pp.normalize_total` followed by
`sc.pp.log1p`.

Each cell is scaled so its total count equals a target. By default the target
is the **median of every cell's total count** (computed over all cells,
including cells with zero counts), reproducing scanpy's `target_sum=None`.
Then `ln(1 + x)` is applied to the scaled values. The sparsity pattern is
preserved exactly because `log1p(0) = 0`.

## Usage

```bash
# scanpy default: median target, log1p
rsomics-sc-normalize filtered_feature_bc_matrix/ -o norm.mtx

# CPM (target 1e6), no log
rsomics-sc-normalize mtx_dir/ --target-sum 1e6 --no-log -o cpm.mtx
```

Input is a 10x MTX directory (`matrix.mtx` or `matrix.mtx.gz`, genes × cells).
Output is a MatrixMarket coordinate matrix in the same genes × cells layout
with real values.

`--exclude-highly-expressed` (scanpy's optional flag, which drops genes
exceeding `max_fraction` of any cell's counts before computing the size
factor) is **not implemented**: it is rarely used in routine pipelines and
adds a two-pass gene-masking step. The size factor here is always the cell's
full total count.

## Origin

This crate is an independent Rust reimplementation of scanpy's
`normalize_total` + `log1p` based on:

- The published method (Wolf, Angerer & Theis, "SCANPY: large-scale
  single-cell gene expression data analysis", *Genome Biology* 2018,
  doi:10.1186/s13059-017-1382-0).
- The public MatrixMarket and 10x Genomics matrix file-format specs.
- Black-box behavior testing against the scanpy Python package.

License: MIT OR Apache-2.0.
Upstream credit: scanpy <https://github.com/scverse/scanpy> (BSD-3-Clause).

#!/usr/bin/env python3
"""scanpy oracle for rsomics-sc-normalize.

Reads a 10x MTX directory, runs sc.pp.normalize_total (default target_sum=None
-> median of per-cell totals) then sc.pp.log1p, and dumps the resulting matrix
as a genes x cells MatrixMarket coordinate file matching the tool's output
layout: header, then `gene cell value` lines in the input's stored order.

Usage: scanpy_normalize_oracle.py <mtx_dir> <out.mtx> [target_sum|median] [--no-log]
"""
import sys

import numpy as np
import scanpy as sc
import scipy.sparse as sp


def main() -> None:
    mtx_dir = sys.argv[1]
    out_path = sys.argv[2]
    target = None
    no_log = "--no-log" in sys.argv[3:]
    for a in sys.argv[3:]:
        if a != "--no-log" and a.lower() != "median":
            target = float(a)

    adata = sc.read_10x_mtx(mtx_dir)
    sc.pp.normalize_total(adata, target_sum=target)
    if not no_log:
        sc.pp.log1p(adata)

    # adata.X is cells x genes; transpose to genes x cells (the input layout).
    gc = sp.csr_matrix(adata.X).T.tocoo()
    n_genes, n_cells = gc.shape
    nnz = gc.nnz

    # Sort to (gene, cell) so the comparison is order-independent of how
    # scipy/scanpy reorders entries internally.
    order = np.lexsort((gc.col, gc.row))
    with open(out_path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{n_genes} {n_cells} {nnz}\n")
        rows = gc.row[order]
        cols = gc.col[order]
        vals = gc.data[order]
        for r, c, v in zip(rows, cols, vals):
            f.write(f"{r + 1} {c + 1} {float(v)!r}\n")


if __name__ == "__main__":
    main()

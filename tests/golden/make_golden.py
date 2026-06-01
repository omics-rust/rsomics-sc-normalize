#!/usr/bin/env python3
"""Generate the deterministic golden 10x MTX dir under tests/golden/tenx.

Genes x cells layout (10x convention). Includes one zero-count cell so the
median definition (over all cells, including zeros) is exercised.
"""
import gzip
import os

import numpy as np
import scipy.io as sio
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "tenx")
os.makedirs(OUT, exist_ok=True)

# Deterministic counts: 60 cells x 40 genes, ~12% density, plus a zero cell.
rng = np.random.default_rng(20260601)
n_cells, n_genes = 60, 40
dense = (rng.random((n_cells, n_genes)) < 0.12) * rng.integers(1, 80, (n_cells, n_genes))
dense[7, :] = 0  # a zero-count cell
Xgc = sp.csr_matrix(dense.T.astype(np.int64))  # genes x cells

mpath = os.path.join(OUT, "matrix.mtx")
sio.mmwrite(mpath, Xgc, field="integer")
with open(mpath, "rb") as f, gzip.open(mpath + ".gz", "wb") as g:
    g.write(f.read())
os.remove(mpath)

with gzip.open(os.path.join(OUT, "barcodes.tsv.gz"), "wt") as f:
    for i in range(n_cells):
        f.write(f"CELL{i:04d}-1\n")
with gzip.open(os.path.join(OUT, "features.tsv.gz"), "wt") as f:
    for i in range(n_genes):
        f.write(f"ENSG{i:08d}\tGene{i}\tGene Expression\n")

print(f"wrote {OUT}: {n_cells} cells x {n_genes} genes, {Xgc.nnz} nonzeros")

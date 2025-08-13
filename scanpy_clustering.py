#!/usr/bin/env python3
"""
General scRNA-seq clustering pipeline using Scanpy.
Performs: PCA → neighbors → Leiden clustering → UMAP → plot.
Usage:
    python scrna_clustering.py data.h5ad output_dir
"""

import sys
import scanpy as sc
import os

adata_file, out_dir = sys.argv[1:]

os.makedirs(out_dir, exist_ok=True)

# === 1. Load data ===
print("Loading data...")
adata = sc.read_h5ad(adata_file)

# === 2. PCA ===
print("PCA...")
sc.tl.pca(adata, svd_solver='arpack')

# === 3. Compute neighbors ===
print("Nearest neighbors...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# === 4. Leiden clustering ===
print("Leiden clustering...")
sc.tl.leiden(adata, resolution=1.0)

# === 5. UMAP embedding ===
print("UMAP...")
sc.tl.umap(adata)

# === 6. Plot UMAP ===
print("Plotting UMAP...")
sc.pl.umap(adata, color=['leiden'], save="_clusters.png", show=False)

# === 7. Save processed AnnData ===
adata.write(os.path.join(out_dir, "adata_clustered.h5ad"))

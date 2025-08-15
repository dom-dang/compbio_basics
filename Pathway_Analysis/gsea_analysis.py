#!/usr/bin/env python3
"""
Run preranked GSEA from DESeq2 results.
Usage:
    python gsea_runner.py deseq_results.csv gene_set_db.gmt out_dir species
Example:
    python gsea_runner.py deseq_results.csv msigdb.gmt gsea_out hsapiens
"""

import sys
import pandas as pd
import gseapy as gp
import os

results_file, gmt_file, out_dir, species = sys.argv[1:]

# === 1. Load DESeq2 results ===
df = pd.read_csv(results_file, index_col=0)

# Require log2FoldChange and p-values
if not {"log2FoldChange", "stat"}.issubset(df.columns):
    raise ValueError("DESeq2 results must contain 'log2FoldChange' and 'stat' columns.")

# === 2. Prepare preranked list ===
# Rank by DESeq2 Wald statistic (stat) or log2FoldChange
rank_series = df["stat"]
rank_series.index.name = "gene"
rank_series = rank_series.sort_values(ascending=False)

# Drop NAs
rank_series = rank_series.dropna()

# === 3. Run GSEA preranked ===
os.makedirs(out_dir, exist_ok=True)

prerank_res = gp.prerank(
    rnk=rank_series,
    gene_sets=gmt_file,  # e.g., MSigDB GMT file
    outdir=out_dir,
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42,
    verbose=True
)

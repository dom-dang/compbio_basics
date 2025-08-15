#!/usr/bin/env python3
"""
Basic QC and filtering for bulk RNA-seq count matrices

Inputs
------
- A raw count matrix (CSV/TSV) with genes as rows, samples as columns.
  The first column must be gene identifiers. Delimiter inferred by file suffix.

Features
--------
- Filter lowly expressed genes by a minimum count in a minimum number of samples.
- Filter samples by a minimum library size (total counts).
- Writes a QC report (YAML) and filtered matrix (CSV).

Examples
--------
python qc_filtering.py counts.csv \
  --min-count 10 --min-samples 3 \
  --min-libsize 1e5 \
  --mito-prefix MT- --max-mito-pct 30 \
  --out-prefix data/filtered

"""
from __future__ import annotations
import argparse, pathlib, sys, math, json
from typing import Tuple, Optional, Iterable
import pandas as pd
import numpy as np
import yaml

def read_table(path: str) -> pd.DataFrame:
    p = pathlib.Path(path)
    if p.suffix.lower() in [".tsv", ".txt"]:
        df = pd.read_csv(p, sep="\t")
    else:
        df = pd.read_csv(p)
    if df.shape[1] < 2:
        raise ValueError("Input must have ≥2 columns: gene_id + ≥1 sample columns.")
    df = df.set_index(df.columns[0])
    return df

def infer_mito_mask(genes: Iterable[str], mito_prefix: Optional[str]) -> np.ndarray:
    if not mito_prefix:
        return np.zeros(len(genes), dtype=bool)
    return pd.Index(genes).str.startswith(mito_prefix).values

def filter_genes(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.Series:
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    return mask

def libsize(series: pd.Series) -> int:
    return int(series.sum())

def mito_pct(counts: pd.DataFrame, mito_mask: np.ndarray) -> pd.Series:
    if mito_mask.sum() == 0:
        return pd.Series(0.0, index=counts.columns, dtype=float)
    mito_counts = counts[mito_mask].sum(axis=0)
    total = counts.sum(axis=0).replace(0, np.nan)
    return (mito_counts / total * 100).fillna(0.0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts", help="Raw count matrix (genes x samples). CSV or TSV.")
    ap.add_argument("--min-count", type=int, default=10, help="Min count per gene to be considered expressed.")
    ap.add_argument("--min-samples", type=int, default=3, help="Min #samples meeting --min-count.")
    ap.add_argument("--min-libsize", type=float, default=0, help="Discard samples with total counts below this.")
    ap.add_argument("--mito-prefix", type=str, default=None, help='Gene prefix for mitochondrial genes (e.g., "MT-").')
    ap.add_argument("--max-mito-pct", type=float, default=100.0, help="Max allowed mitochondrial %% per sample.")
    ap.add_argument("--out-prefix", required=True, help="Output path prefix (no extension).")
    args = ap.parse_args()

    counts = read_table(args.counts).astype(np.int64)
    mito_mask = infer_mito_mask(counts.index, args.mito_prefix)

    # Sample QC
    libsizes = counts.sum(axis=0)
    mito_perc = mito_pct(counts, mito_mask)
    keep_samples = (libsizes >= args.min_libsize) & (mito_perc <= args.max_mito_pct)
    dropped_samples = counts.columns[~keep_samples].tolist()
    counts = counts.loc[:, keep_samples]

    # Gene QC
    keep_genes = filter_genes(counts, args.min_count, args.min_samples)
    dropped_genes = counts.index[~keep_genes].tolist()
    counts = counts.loc[keep_genes]

    # Write outputs
    out_mat = f"{args.out_prefix}.filtered_counts.csv"
    counts.to_csv(out_mat)

    qc = {
        "input": pathlib.Path(args.counts).as_posix(),
        "output_matrix": out_mat,
        "n_samples_in": int(libsizes.shape[0]),
        "n_samples_out": int(keep_samples.sum()),
        "dropped_samples": dropped_samples,
        "min_libsize": args.min_libsize,
        "mito_prefix": args.mito_prefix,
        "max_mito_pct": args.max_mito_pct,
        "sample_libsizes": libsizes.to_dict(),
        "sample_mito_pct": mito_perc.round(3).to_dict(),
        "n_genes_in": int(keep_genes.shape[0] + len(dropped_genes)),
        "n_genes_out": int(keep_genes.sum()),
        "dropped_genes": dropped_genes[:50],  # preview
        "gene_filter_min_count": args.min_count,
        "gene_filter_min_samples": args.min_samples,
    }
    with open(f"{args.out_prefix}.qc.yaml", "w") as f:
        yaml.safe_dump(qc, f, sort_keys=False)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Normalize bulk RNA-seq count matrices

Inputs
------
- Count matrix (CSV/TSV) with genes as rows, samples as columns.
- Optional gene length file (CSV/TSV) with columns: gene_id, length_bp (for TPM).

Methods
-------
- CPM: counts per million.
- TPM: requires gene lengths (bp or kb); internally converts to kb.
- MOR: Median-of-ratios (DESeq-style size factors) on geometric means.

Options
-------
--method {cpm,tpm,mor}
--log1p  (apply log(1+x) after normalization)
--scale  (z-score genes across samples)

Examples
--------
python normalization.py data/filtered_counts.csv --method cpm --log1p \
  --out data/filtered.cpm.log1p.csv

python normalization.py data/filtered_counts.csv --method tpm \
  --gene-lengths gene_lengths.tsv --out data/filtered.tpm.csv
"""
from __future__ import annotations
import argparse, pathlib
import numpy as np
import pandas as pd

def read_table(path: str) -> pd.DataFrame:
    p = pathlib.Path(path)
    sep = "\t" if p.suffix.lower() in [".tsv", ".txt"] else ","
    df = pd.read_csv(p, sep=sep).set_index(lambda c: c if c != 0 else None)
    if df.index.name is None:
        df = df.set_index(df.columns[0])
    return df

def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    libsize = counts.sum(axis=0)
    return counts.divide(libsize, axis=1) * 1e6

def tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    lengths_kb = lengths_bp / 1_000.0
    # Align
    common = counts.index.intersection(lengths_kb.index)
    c = counts.loc[common]
    lkb = lengths_kb.loc[common].replace(0, np.nan)
    rpk = c.divide(lkb, axis=0)
    per_sample_scaler = rpk.sum(axis=0) / 1e6
    return rpk.divide(per_sample_scaler, axis=1)

def mor(counts: pd.DataFrame) -> pd.DataFrame:
    # Geometric means across samples (ignore zeros via log trick)
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts.replace(0, np.nan))
    geo_means = np.exp(np.nanmean(log_counts, axis=1))
    gm = pd.Series(geo_means, index=counts.index).replace(0, np.nan)
    ratios = counts.divide(gm, axis=0).replace([np.inf, -np.inf], np.nan)
    size_factors = ratios.median(axis=0, skipna=True)
    norm = counts.divide(size_factors, axis=1)
    return norm

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts", help="Counts matrix (genes x samples).")
    ap.add_argument("--method", choices=["cpm","tpm","mor"], default="cpm")
    ap.add_argument("--gene-lengths", help="Gene length table (gene_id, length_bp) for TPM.")
    ap.add_argument("--log1p", action="store_true", help="Apply log1p after normalization.")
    ap.add_argument("--scale", action="store_true", help="Z-score genes across samples after other steps.")
    ap.add_argument("--out", required=True, help="Output CSV path.")
    args = ap.parse_args()

    counts = read_table(args.counts).astype(float)

    if args.method == "cpm":
        mat = cpm(counts)
    elif args.method == "tpm":
        if not args.gene_lengths:
            raise SystemExit("--gene-lengths is required for TPM.")
        gl = read_table(args.gene_lengths)
        # Accept any column name for length
        length_col = gl.columns[0]
        mat = tpm(counts, gl[length_col].astype(float))
    else:
        mat = mor(counts)

    if args.log1p:
        mat = np.log1p(mat)

    if args.scale:
        mat = (mat - mat.mean(axis=1).values.reshape(-1,1)) / mat.std(axis=1).replace(0, np.nan).values.reshape(-1,1)

    mat.to_csv(args.out)

if __name__ == "__main__":
    main()

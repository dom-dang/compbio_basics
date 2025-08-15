#!/usr/bin/env python3
"""
Plots for GO/GSEA results you already generate

Assumes a tabular enrichment file with columns like:
- term / pathway / description
- pval or pvalue (optional)
- padj or qvalue or fdr
- NES or enrichment_score (optional)
- gene_ratio (e.g., "12/250") or overlap (e.g., "12/250")

Creates barplot of top N by -log10(adjusted p) and optional dotplot.

Examples
--------
python pathway_visualization.py go_results.csv --top 15 --bar go_top15.png
python pathway_visualization.py gsea_results.tsv --top 20 --dot gsea_dot.png
"""
import argparse, math, pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_table(path: str) -> pd.DataFrame:
    p = pathlib.Path(path)
    sep = "\t" if p.suffix.lower() in [".tsv", ".txt"] else ","
    return pd.read_csv(p, sep=sep)

def coalesce(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            return df[c]
    raise KeyError(f"None of {cols} found in columns: {df.columns.tolist()}")

def parse_ratio(s: str) -> float:
    if pd.isna(s): return np.nan
    if isinstance(s, (int, float)): return float(s)
    s = str(s)
    if "/" in s:
        a, b = s.split("/", 1)
        a = float(a); b = float(b) if float(b) != 0 else np.nan
        return a / b
    try:
        return float(s)
    except Exception:
        return np.nan

def make_bar(df: pd.DataFrame, out_path: str, top: int):
    df = df.sort_values("-log10padj", ascending=False).head(top).iloc[::-1]
    plt.figure(figsize=(8, max(4, 0.4 * len(df))))
    plt.barh(df["term"], df["-log10padj"])
    plt.xlabel("-log10(adjusted p)")
    plt.ylabel("Pathway")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def make_dot(df: pd.DataFrame, out_path: str, top: int):
    df = df.sort_values("-log10padj", ascending=False).head(top).iloc[::-1]
    sizes = 200 * df["gene_ratio"].fillna(0.0)
    colors = df.get("NES", df.get("enrichment_score", pd.Series(np.nan, index=df.index)))
    plt.figure(figsize=(8, max(4, 0.45 * len(df))))
    sc = plt.scatter(df["-log10padj"], df["term"], s=sizes, c=colors)
    plt.xlabel("-log10(adjusted p)")
    plt.ylabel("Pathway")
    if not colors.isna().all():
        cbar = plt.colorbar(sc)
        cbar.set_label("NES / enrichment score")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("enrichment", help="GO/GSEA results table (CSV/TSV).")
    ap.add_argument("--top", type=int, default=20, help="Top N pathways to visualize.")
    ap.add_argument("--bar", help="Output PNG for horizontal bar plot.")
    ap.add_argument("--dot", help="Output PNG for dot plot (size ~ gene ratio; color ~ NES/ES if available).")
    args = ap.parse_args()

    df = read_table(args.enrichment)
    # Standardize key fields
    df = df.rename(columns={
        "description":"term", "pathway":"term", "term_name":"term", "ID":"term"
    })
    if "term" not in df.columns:
        # fall back to first column as label
        df = df.rename(columns={df.columns[0]:"term"})

    padj = coalesce(df, ["padj","qvalue","fdr","adj_p","adj_pval","p_adjusted"])
    df["-log10padj"] = -np.log10(padj.astype(float).replace(0, np.nextafter(0, 1)))
    if "gene_ratio" in df.columns:
        df["gene_ratio"] = df["gene_ratio"].map(parse_ratio)
    elif "overlap" in df.columns:
        df["gene_ratio"] = df["overlap"].map(parse_ratio)
    else:
        df["gene_ratio"] = np.nan

    if "NES" not in df.columns and "enrichment_score" not in df.columns:
        # ok, color will be absent
        pass

    if args.bar:
        make_bar(df.copy(), args.bar, args.top)
    if args.dot:
        make_dot(df.copy(), args.dot, args.top)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Compare gene lists, compute pairwise overlaps, Jaccard index,
and Fisher's exact test p-values. Optionally produce a 2–3-set Venn diagram.

Inputs
------
- One or more plain-text files, each containing a gene per line.
- Either provide a background/universe size with --background N, or a file with
  the universe genes using --universe file.txt (recommended).

Outputs
-------
- CSV matrices: overlap_counts.csv, jaccard.csv, fisher_p.csv
- Optional PNG Venn (for 2 or 3 lists).

Example
-------
python gene_set_overlap.py a.txt b.txt c.txt --universe all_genes.txt --venn venn.png
"""
import argparse, itertools, math, pathlib
from typing import List, Dict, Set
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt

def read_list(path: str) -> Set[str]:
    with open(path) as f:
        return set(x.strip() for x in f if x.strip())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lists", nargs="+", help="Gene list files (txt, one gene per line).")
    ap.add_argument("--universe", help="File containing background gene universe (one per line).")
    ap.add_argument("--background", type=int, help="Size of background, if --universe not provided.")
    ap.add_argument("--venn", help="Output PNG path for Venn diagram (only if 2 or 3 lists).")
    ap.add_argument("--out-prefix", default="overlap", help="Prefix for CSV outputs.")
    args = ap.parse_args()

    names = [pathlib.Path(p).stem for p in args.lists]
    sets = [read_list(p) for p in args.lists]

    if args.universe:
        universe = read_list(args.universe)
        bg_size = len(universe)
        sets = [s & universe for s in sets]
    elif args.background:
        bg_size = args.background
    else:
        # fallback to union as background
        universe = set().union(*sets)
        bg_size = len(universe)

    n = len(sets)
    overlap = pd.DataFrame(0, index=names, columns=names, dtype=int)
    jacc = pd.DataFrame(0.0, index=names, columns=names, dtype=float)
    pvals = pd.DataFrame(1.0, index=names, columns=names, dtype=float)

    for i in range(n):
        for j in range(n):
            inter = len(sets[i] & sets[j])
            union = len(sets[i] | sets[j])
            overlap.iat[i, j] = inter
            jacc.iat[i, j] = inter / union if union else 0.0
            if i != j:
                a = inter
                b = len(sets[i]) - inter
                c = len(sets[j]) - inter
                d = bg_size - (a + b + c)
                _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                pvals.iat[i, j] = p

    overlap.to_csv(f"{args.out-prefix}.overlap_counts.csv")
    jacc.to_csv(f"{args.out-prefix}.jaccard.csv")
    pvals.to_csv(f"{args.out-prefix}.fisher_p.csv")

    if args.venn and n in (2,3):
        from matplotlib_venn import venn2, venn3  # requires matplotlib-venn
        plt.figure()
        if n == 2:
            venn2(subsets=(sets[0], sets[1]), set_labels=(names[0], names[1]))
        else:
            venn3(subsets=(sets[0], sets[1], sets[2]), set_labels=(names[0], names[1], names[2]))
        plt.tight_layout()
        plt.savefig(args.venn, dpi=200)
        plt.close()

if __name__ == "__main__":
    main()

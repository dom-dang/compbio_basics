#!/usr/bin/env python3
"""
Make a volcano plot from DESeq2 results
Usage:
    python volcano_plot.py deseq_results.csv
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

results_file = sys.argv[1]
df = pd.read_csv(results_file, index_col=0)

plt.figure(figsize=(7,6))
plt.scatter(df['log2FoldChange'], -np.log10(df['padj']), 
            c=(df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1),
            cmap="coolwarm", alpha=0.7)

plt.axhline(-np.log10(0.05), color="grey", linestyle="--")
plt.axvline(1, color="grey", linestyle="--")
plt.axvline(-1, color="grey", linestyle="--")

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 Adjusted P-value")
plt.tight_layout()

# plt.savefig("volcano_plot.png", dpi=300)


#!/usr/bin/env python3
"""
Plot PCA from normalized counts
Usage:
    python pca_analysis.py normalized_counts.csv metadata.csv condition
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

counts_file, metadata_file, condition_col = sys.argv[1:]

counts_df = pd.read_csv(counts_file, index_col=0)
metadata_df = pd.read_csv(metadata_file, index_col=0).loc[counts_df.columns]

pca = PCA(n_components=2)
pcs = pca.fit_transform(counts_df.T)

plt.figure(figsize=(6,6))
for cond in metadata_df[condition_col].unique():
    idx = metadata_df[condition_col] == cond
    plt.scatter(pcs[idx,0], pcs[idx,1], label=cond)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
plt.legend()
plt.tight_layout()

# plt.savefig("pca_plot.png", dpi=300)

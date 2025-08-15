#!/usr/bin/env python3
"""
GO enrichment from DESeq2 significant genes
Usage:
    python go_enrichment.py deseq_results.csv species_code
Example:
    python go_enrichment.py deseq_results.csv hsapiens
"""

import sys
import pandas as pd
from gprofiler import GProfiler

results_file, species_code = sys.argv[1:]
df = pd.read_csv(results_file, index_col=0)

sig_genes = df.query("padj < 0.05").index.tolist()

gp = GProfiler(return_dataframe=True)
res = gp.profile(organism=species_code, query=sig_genes)

# res.to_csv("go_enrichment.csv", index=False)


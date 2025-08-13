#!/usr/bin/env python3
"""
General DESeq2 pipeline using PyDESeq2
Usage:
    python deseq2_pipeline.py counts.csv metadata.csv condition control treatment
"""

import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts_file, metadata_file, condition_col, ref_level, alt_level = sys.argv[1:]

counts_df = pd.read_csv(counts_file, index_col=0)
metadata_df = pd.read_csv(metadata_file, index_col=0).loc[counts_df.columns]

dds = DeseqDataSet(
    counts=counts_df,
    clinical=metadata_df,
    design_factors=[condition_col],
    ref_level={condition_col: ref_level}
)

dds.deseq2()

stat_res = DeseqStats(dds, contrast=[condition_col, alt_level, ref_level])
stat_res.summary()

# stat_res.results_df.to_csv("deseq_results.csv")
# dds.norm_counts.to_csv("normalized_counts.csv")



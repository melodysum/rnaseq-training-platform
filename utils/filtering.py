"""
utils/filtering.py
------------------
Low-expression filtering logic and related statistics.
"""

import numpy as np
import pandas as pd


def filter_low_expression(counts: pd.DataFrame,
                           min_count: int = 10,
                           min_samples: int = 5) -> pd.DataFrame:
    """Keep genes with >= min_count in at least min_samples samples."""
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[mask]


def expression_summary(counts: pd.DataFrame) -> pd.DataFrame:
    """Per-gene summary: mean, median, max count, fraction of zero samples."""
    return pd.DataFrame({
        "mean_count":    counts.mean(axis=1),
        "median_count":  counts.median(axis=1),
        "max_count":     counts.max(axis=1),
        "pct_zeros":     (counts == 0).mean(axis=1) * 100,
    })


def check_gene_fate(counts_raw: pd.DataFrame,
                    counts_filtered: pd.DataFrame,
                    gene_list: list) -> pd.DataFrame:
    """
    For each gene in gene_list, report whether it was:
      - retained after filtering
      - removed due to low expression
      - not found in the dataset
    """
    rows = []
    all_genes      = set(counts_raw.index)
    retained_genes = set(counts_filtered.index)

    for gene in gene_list:
        gene = gene.strip()
        if gene not in all_genes:
            status = "❓ Not found in dataset"
            mean_count = None
        elif gene in retained_genes:
            status = "✅ Retained after filtering"
            mean_count = counts_raw.loc[gene].mean()
        else:
            status = "🔴 Removed (low expression)"
            mean_count = counts_raw.loc[gene].mean()

        rows.append({"Gene": gene, "Status": status,
                     "Mean count (all samples)": mean_count})

    return pd.DataFrame(rows)


def normalise_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Convert raw counts to CPM."""
    lib_sizes = counts.sum(axis=0)
    return counts.divide(lib_sizes, axis=1) * 1e6


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log2(CPM + 1)."""
    return np.log2(normalise_cpm(counts) + 1)

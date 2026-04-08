"""
utils/data_loader.py
--------------------
Handles loading, validation, and demo-data fallback for the RNA-seq app.
"""

import os
import numpy as np
import pandas as pd
import streamlit as st

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def load_demo_data():
    """Load the built-in demo dataset (real counts + metadata)."""
    counts = pd.read_csv(os.path.join(DATA_DIR, "counts.csv"), index_col=0)
    metadata = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"), index_col=0)
    return counts, metadata


def validate_and_parse(counts_file, metadata_file):
    """
    Parse and validate uploaded CSV files with comprehensive checks.
    Returns (counts_df, metadata_df, error_message).
    error_message is None if everything passed.

    Input expectations:
    - counts.csv : genes as rows, samples as columns, raw integer counts
    - metadata.csv: samples as rows, must include 'groupA' column
    - Input must be raw counts, NOT TPM / FPKM / CPM
    """
    # ── 1. Read CSVs ──────────────────────────────────────────────────────────
    try:
        counts = pd.read_csv(counts_file, index_col=0)
    except Exception as e:
        return None, None, f"Could not read counts.csv: {e}"

    try:
        metadata = pd.read_csv(metadata_file, index_col=0)
    except Exception as e:
        return None, None, f"Could not read metadata.csv: {e}"

    # ── 2. Empty / all-NA check ───────────────────────────────────────────────
    if counts.empty or counts.isna().all().all():
        return None, None, (
            "counts.csv appears to be empty or contains only missing values."
        )
    if metadata.empty or metadata.isna().all().all():
        return None, None, (
            "metadata.csv appears to be empty or contains only missing values."
        )

    # ── 3. Minimum dimensions ─────────────────────────────────────────────────
    if counts.shape[1] < 2:
        return None, None, (
            f"counts.csv must have at least 2 sample columns. "
            f"Found {counts.shape[1]}. Check that gene names are in the first "
            f"column (used as the row index) and each remaining column is a sample."
        )
    if counts.shape[0] < 1:
        return None, None, "counts.csv must contain at least 1 gene row."
    if len(metadata) < 2:
        return None, None, (
            f"metadata.csv must have at least 2 samples. Found {len(metadata)}."
        )

    # ── 4. Required metadata column ───────────────────────────────────────────
    if "groupA" not in metadata.columns:
        return None, None, (
            "metadata.csv must contain a 'groupA' column identifying sample groups. "
            "Optional but important: 'donor' (paired DE), 'batch' (batch correction), "
            "'sex', 'age' (descriptive)."
        )

    # ── 5. Duplicate checks ───────────────────────────────────────────────────
    dup_count_cols = counts.columns[counts.columns.duplicated()].tolist()
    if dup_count_cols:
        return None, None, (
            f"Duplicate sample names found in counts.csv columns: "
            f"{', '.join(dup_count_cols)}. Each sample must have a unique name."
        )
    dup_meta_idx = metadata.index[metadata.index.duplicated()].tolist()
    if dup_meta_idx:
        return None, None, (
            f"Duplicate sample names found in metadata.csv index: "
            f"{', '.join(str(d) for d in dup_meta_idx)}."
        )
    dup_genes = counts.index[counts.index.duplicated()].tolist()
    if dup_genes:
        preview = ', '.join(str(g) for g in dup_genes[:5])
        suffix = '...' if len(dup_genes) > 5 else ''
        return None, None, (
            f"Duplicate gene identifiers found in counts.csv index: "
            f"{preview}{suffix}. Each gene must appear only once."
        )

    # ── 6. Sample name alignment ──────────────────────────────────────────────
    count_samples    = set(counts.columns)
    meta_samples     = set(metadata.index)
    missing_in_meta  = count_samples - meta_samples
    missing_in_count = meta_samples  - count_samples

    if missing_in_meta:
        return None, None, (
            f"These sample names in counts.csv are missing from metadata.csv: "
            f"{', '.join(sorted(missing_in_meta))}. "
            "Sample names must match exactly (case-sensitive)."
        )
    if missing_in_count:
        return None, None, (
            f"These samples in metadata.csv are not found in counts.csv columns: "
            f"{', '.join(sorted(missing_in_count))}. "
            "All metadata samples must be present in the count matrix."
        )

    # ── 7. Numeric counts check ───────────────────────────────────────────────
    try:
        counts = counts.apply(pd.to_numeric, errors="raise")
    except Exception:
        return None, None, (
            "counts.csv contains non-numeric values. "
            "All count values must be integers or floats. "
            "Ensure the first column contains gene identifiers (row index) "
            "and all other columns contain numeric counts."
        )

    # ── 8. Negative values ────────────────────────────────────────────────────
    if (counts < 0).any().any():
        return None, None, (
            "counts.csv contains negative values. "
            "Raw counts must be non-negative integers. "
            "This app requires raw counts — not TPM, FPKM, or CPM values, "
            "which can produce negative log-transformed values."
        )

    # ── 9. Reorder metadata to match counts column order ─────────────────────
    metadata = metadata.reindex(counts.columns)

    return counts, metadata, None


def get_paired_columns(metadata: pd.DataFrame):
    """Return (donors, ctrl_cols, treat_cols) for paired design."""
    if "donor" not in metadata.columns:
        return None, None, None
    donors = metadata["donor"].unique()
    ctrl_cols  = [f"{d}_control"   for d in donors if f"{d}_control"   in metadata.index]
    treat_cols = [f"{d}_treatment" for d in donors if f"{d}_treatment" in metadata.index]
    return donors, ctrl_cols, treat_cols


def init_session_data():
    """
    Initialise session_state with demo data if nothing uploaded yet.
    Call once at the top of each page.
    """
    if "counts" not in st.session_state:
        counts, metadata = load_demo_data()
        st.session_state["counts"]      = counts
        st.session_state["metadata"]    = metadata
        st.session_state["data_source"] = "demo"

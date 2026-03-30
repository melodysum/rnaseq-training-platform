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
    Parse and validate uploaded CSV files.
    Returns (counts_df, metadata_df, error_message).
    error_message is None if everything is fine.
    """
    try:
        counts = pd.read_csv(counts_file, index_col=0)
    except Exception as e:
        return None, None, f"Could not read counts.csv: {e}"

    try:
        metadata = pd.read_csv(metadata_file, index_col=0)
    except Exception as e:
        return None, None, f"Could not read metadata.csv: {e}"

    # Check required metadata column
    if "groupA" not in metadata.columns:
        return None, None, "metadata.csv must contain a 'groupA' column."

    # Check sample name alignment
    count_samples = set(counts.columns)
    meta_samples  = set(metadata.index)
    missing_in_meta  = count_samples - meta_samples
    missing_in_count = meta_samples  - count_samples

    if missing_in_meta:
        return None, None, (
            f"These sample columns in counts.csv are not in metadata.csv: "
            f"{', '.join(sorted(missing_in_meta))}"
        )
    if missing_in_count:
        return None, None, (
            f"These samples in metadata.csv are missing from counts.csv: "
            f"{', '.join(sorted(missing_in_count))}"
        )

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
        st.session_state["counts"]   = counts
        st.session_state["metadata"] = metadata
        st.session_state["data_source"] = "demo"

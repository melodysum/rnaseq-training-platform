"""
pages/10_Public_Data.py — Lesson 10: Public Data & Reproducible Analysis (Coming soon)
"""
import streamlit as st
from utils.data_loader import init_session_data

st.set_page_config(page_title="Lesson 10 — Public Data & Reproducibility", page_icon="🌐", layout="wide")
init_session_data()

st.title("🌐 Lesson 10 — Public Data & Reproducible Analysis")
st.markdown("> This lesson is currently being developed.")
st.info("""
**Coming soon — what this lesson will cover:**
- Downloading RNA-seq datasets from GEO and ArrayExpress
- Understanding GEO metadata and supplementary files
- Ensuring reproducible workflows with version control
- Conda environments, Snakemake, and Nextflow concepts
- Sharing and documenting analysis for publication
""")
st.divider()
st.page_link("pages/9_Functional_Enrichment.py", label="← Lesson 9: Functional Enrichment", icon="🧩")

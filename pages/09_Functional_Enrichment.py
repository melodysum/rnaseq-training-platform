"""
pages/09_Functional_Enrichment.py — Lesson 9: Functional Enrichment (Coming soon)
"""
import streamlit as st
from utils.data_loader import init_session_data

st.set_page_config(page_title="Lesson 9 — Functional Enrichment", page_icon="🧩", layout="wide")
init_session_data()

st.title("🧩 Lesson 9 — Functional Enrichment (GO / GSEA)")
st.markdown("> This lesson is currently being developed.")
st.info("""
**Coming soon — what this lesson will cover:**
- Gene Ontology (GO) term analysis
- Over-representation analysis (ORA)
- Gene Set Enrichment Analysis (GSEA)
- Pathway databases: KEGG, Reactome, MSigDB
- Interpreting and visualising enrichment results
""")
st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/08_Clustering_Heatmaps.py", label="← Lesson 8: Clustering & Heatmaps", icon="🗺️")
with col_n2:
    st.page_link("pages/10_Public_Data.py", label="Lesson 10: Public Data →", icon="🌐")

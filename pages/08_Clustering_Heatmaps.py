"""
pages/08_Clustering_Heatmaps.py — Lesson 8: Clustering & Heatmaps (Coming soon)
"""
import streamlit as st
from utils.data_loader import init_session_data

st.set_page_config(page_title="Lesson 8 — Clustering & Heatmaps", page_icon="🗺️", layout="wide")
init_session_data()

st.title("🗺️ Lesson 8 — Clustering & Heatmaps")
st.markdown("> This lesson is currently being developed.")
st.info("""
**Coming soon — what this lesson will cover:**
- Hierarchical clustering of genes and samples
- Heatmap interpretation
- Choosing clustering methods and distance metrics
- Identifying gene modules and expression patterns
- Integration with DE results
""")
st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/07_Differential_Expression.py", label="← Lesson 7: Differential Expression", icon="🧪")
with col_n2:
    st.page_link("pages/09_Functional_Enrichment.py", label="Lesson 9: Functional Enrichment →", icon="🧩")

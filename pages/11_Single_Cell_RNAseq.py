"""
pages/11_Single_Cell_RNAseq.py
Lesson 11 — Single-Cell RNA-seq (Introductory Module)

Conceptual introduction only. No real scRNA-seq pipeline is implemented.
All visualisations use simulated data for illustration purposes.
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(
    page_title="Lesson 11 — Single-Cell RNA-seq",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 Lesson 11 — Single-Cell RNA-seq (Introductory Module)")
st.info(
    "**This is a conceptual introduction only.** "
    "This app does not implement a real single-cell analysis pipeline. "
    "All visualisations use simulated data for illustration.",
    icon="ℹ️",
)

# ════════════════════════════════════════════════════════════════════════════════
# 1. BULK vs SINGLE-CELL
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("📊 Bulk RNA-seq vs Single-Cell RNA-seq")

comparison_data = {
    "Property": [
        "Unit of measurement",
        "Biological resolution",
        "Typical output",
        "Sparsity",
        "Common visualisations",
        "Sample sizes",
        "Computational complexity",
        "Typical use cases",
    ],
    "Bulk RNA-seq": [
        "Average across thousands–millions of cells per sample",
        "Sample-level; obscures cell-to-cell heterogeneity",
        "Gene × sample count matrix (dense)",
        "Low — most genes are detected per sample",
        "PCA, heatmaps, volcano plots, MA plots",
        "5–100+ samples (individuals, time points, conditions)",
        "Moderate — standard R/Python pipelines",
        "Tissue-level DE, population studies, biomarker discovery",
    ],
    "Single-Cell RNA-seq": [
        "Individual cell transcriptome",
        "Cell-level; reveals subpopulations and rare cell types",
        "Gene × cell count matrix (very sparse)",
        "High — most genes have zero counts per cell (dropouts)",
        "UMAP, t-SNE, violin plots, dot plots, trajectory plots",
        "Thousands to millions of cells; fewer samples",
        "High — large matrices, graph-based clustering, integration",
        "Cell-type discovery, trajectory analysis, cell–cell communication",
    ],
}

st.dataframe(
    pd.DataFrame(comparison_data).set_index("Property"),
    use_container_width=True,
)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 2. KEY CONCEPTS
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("💡 Key Single-Cell Concepts")

col1, col2 = st.columns(2)
with col1:
    st.markdown("""
### Cells vs Samples
In bulk RNA-seq, each **column** of your count matrix represents a biological
**sample** (e.g. a tissue biopsy from one patient). The signal you measure is
the average expression across all cells in that tissue.

In scRNA-seq, each **column** represents a single **cell**. A typical 10x
Chromium experiment captures 2,000–10,000 cells. The result is a huge matrix
— often 20,000 genes × 5,000 cells — where most entries are zero.

### Sparsity and Dropouts
Sequencing depth per cell is much lower than in bulk RNA-seq. Many transcripts
are simply not captured even if they are present in the cell. These **dropouts**
result in a highly sparse matrix.

This sparsity is not purely noise — it reflects the stochastic nature of mRNA
capture. Special normalisation methods (e.g. scran, SCTransform) are designed
to handle this, rather than simply applying CPM or log-normalisation.
    """)
with col2:
    st.markdown("""
### Clustering and Cell Types
Because cells within a tissue are heterogeneous, the first goal of scRNA-seq
is usually to **group cells by transcriptional similarity**. Graph-based
clustering (Leiden, Louvain) on a shared-nearest-neighbour graph is standard.

Each cluster (ideally) corresponds to a cell population — but clusters are
not inherently labelled. You identify what each cluster is by examining its
**marker genes**: genes highly and specifically expressed in that cluster
compared to all others.

### UMAP
UMAP (Uniform Manifold Approximation and Projection) compresses the
high-dimensional transcriptomic space (thousands of genes) into 2D for
visualisation. Points close together in UMAP space are transcriptionally
similar. UMAP is used to visualise cluster structure, trajectory, and
marker gene expression.

**Important:** UMAP coordinates are not interpretable as distances — only
the neighbourhood relationships matter.

### Batch Effects and Integration
When cells from multiple experiments, donors, or sequencing runs are combined,
technical variation (batch effects) can dominate over biology. Integration
methods (Harmony, Seurat CCA, scVI) align datasets so that cell types from
different batches cluster together rather than by batch.
    """)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 3. SIMULATED UMAP VISUALISATION
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🗺️ Conceptual UMAP Visualisation")
st.caption("**Simulated data for illustration only.** Not real scRNA-seq data.")

# Reproducible simulated clusters
rng = np.random.default_rng(42)

CLUSTER_DEFS = [
    ("T cells",          (-4.0, 3.5),   350,  "#2563eb"),
    ("Monocytes",        (4.5, 2.0),    280,  "#dc2626"),
    ("NK cells",         (-3.5, -3.0),  180,  "#16a34a"),
    ("B cells",          (1.5, -4.5),   220,  "#9333ea"),
    ("Epithelial cells", (6.0, -1.5),   160,  "#ea580c"),
    ("Macrophages",      (2.5, 4.5),    200,  "#0891b2"),
    ("Fibroblasts",      (-1.0, 0.5),    90,  "#854d0e"),
]

rows = []
for name, (cx, cy), n, color in CLUSTER_DEFS:
    xs = rng.normal(cx, 0.9, n)
    ys = rng.normal(cy, 0.9, n)
    for x, y in zip(xs, ys):
        rows.append({"UMAP_1": x, "UMAP_2": y, "Cluster": name})

umap_df = pd.DataFrame(rows)

fig_umap = px.scatter(
    umap_df, x="UMAP_1", y="UMAP_2", color="Cluster",
    title="Simulated UMAP — 7 Toy Cell Clusters<br><sup>Simulated data for illustration only</sup>",
    color_discrete_map={name: color for name, _, _, color in CLUSTER_DEFS},
    opacity=0.7,
    labels={"UMAP_1": "UMAP 1", "UMAP_2": "UMAP 2"},
    height=500,
)
fig_umap.update_traces(marker=dict(size=4))
fig_umap.update_layout(legend=dict(title="Cell Cluster", itemsizing="constant"))
st.plotly_chart(fig_umap, use_container_width=True)

# Cluster marker summary table
st.markdown("#### Example Marker Genes per Cluster")
st.caption("Hardcoded illustrative examples — real marker genes from published literature.")

marker_table = pd.DataFrame([
    {"Cluster": "T cells",          "Key Markers": "CD3D, CD3E, IL7R, LTB, TRAC",
     "Notes": "Pan-T cell markers; IL7R marks naive/memory T cells"},
    {"Cluster": "Monocytes",        "Key Markers": "LST1, S100A8, S100A9, CTSS, FCN1",
     "Notes": "Classical monocyte markers; S100A8/A9 are danger signals"},
    {"Cluster": "NK cells",         "Key Markers": "GNLY, NKG7, GZMB, FCER1G, KLRD1",
     "Notes": "Cytotoxic NK cell markers; GNLY/NKG7 mark cytotoxic granules"},
    {"Cluster": "B cells",          "Key Markers": "CD79A, MS4A1, BANK1, CD22, HLA-DQA1",
     "Notes": "B cell identity; MS4A1 = CD20"},
    {"Cluster": "Epithelial cells", "Key Markers": "EPCAM, KRT8, KRT18, KRT19, CLDN4",
     "Notes": "Epithelial identity; keratins define epithelial origin"},
    {"Cluster": "Macrophages",      "Key Markers": "CD68, MRC1, MARCO, APOE, C1QA",
     "Notes": "Tissue-resident macrophage markers; C1Q marks mature macrophages"},
    {"Cluster": "Fibroblasts",      "Key Markers": "COL1A1, COL1A2, DCN, LUM, THY1",
     "Notes": "Stromal fibroblast markers; collagens define ECM-producing cells"},
])
st.dataframe(marker_table.set_index("Cluster"), use_container_width=True)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 4. TYPICAL scRNA-seq WORKFLOW
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🔧 Typical Single-Cell Analysis Workflow")
st.caption("Conceptual overview — not implemented in this app.")

workflow_steps = [
    ("1. Raw data & alignment",
     "FASTQ → Cell Ranger (10x) or STARsolo → filtered_feature_bc_matrix. "
     "Output: cells × genes UMI count matrix."),
    ("2. Quality control",
     "Filter cells by: min genes detected, max % mitochondrial reads (marks dying cells), "
     "min UMI count. Remove doublets with tools like DoubletFinder or Scrublet."),
    ("3. Normalisation",
     "Normalise by library size (scran pooling, or simple CPM). Log-transform. "
     "SCTransform (Seurat) models technical variation using negative binomial regression."),
    ("4. Feature selection",
     "Select highly variable genes (HVGs) — typically 2,000–5,000 genes — "
     "to reduce noise and focus on informative variation."),
    ("5. Dimensionality reduction",
     "PCA on HVGs (typically 10–50 PCs retained). Then UMAP or t-SNE "
     "on the PCA embedding for 2D visualisation."),
    ("6. Clustering",
     "Build k-nearest-neighbour graph on PCA coordinates. "
     "Apply Leiden or Louvain community detection. Resolution parameter controls granularity."),
    ("7. Marker gene identification",
     "Find genes differentially expressed between each cluster and all others "
     "(Wilcoxon rank-sum test). Top markers define cluster identity."),
    ("8. Cell-type annotation",
     "Match marker genes to known cell-type signatures from literature or "
     "reference datasets (e.g. HCA, CellTypist, SingleR)."),
    ("9. Downstream analysis",
     "Differential abundance, trajectory (pseudotime), cell–cell communication, "
     "multi-sample integration, or integration with spatial data."),
]

for step, desc in workflow_steps:
    with st.expander(step):
        st.markdown(desc)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 5. WHAT THIS APP DOES NOT IMPLEMENT
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🚧 What This App Does Not Yet Implement")
st.error(
    "The following single-cell capabilities are **not available** in this web app. "
    "Use Scanpy (Python) or Seurat (R) for real scRNA-seq analysis.",
    icon="🚧",
)

not_implemented = [
    "Cell-level QC filtering (mitochondrial reads, doublet removal)",
    "Single-cell normalisation (scran, SCTransform, or similar)",
    "Highly variable gene (HVG) selection",
    "PCA and shared nearest-neighbour (SNN) graph on real scRNA-seq data",
    "True UMAP on real single-cell matrices",
    "Leiden/Louvain graph-based clustering on real data",
    "Cluster-level differential expression",
    "Cell-type annotation from marker genes",
    "Multi-sample integration (Harmony, scVI, Seurat CCA)",
    "Trajectory / pseudotime analysis (Monocle, scVelo)",
    "AnnData (.h5ad) or Seurat (.rds) file format support",
    "10x Cell Ranger output parsing",
]

for item in not_implemented:
    st.markdown(f"- ❌ {item}")

st.markdown("""
**Recommended tools for real scRNA-seq analysis:**
- **Python:** [Scanpy](https://scanpy.readthedocs.io/) + AnnData
- **R:** [Seurat](https://satijalab.org/seurat/) or [Bioconductor scran/scater](https://bioconductor.org/books/release/OSCA/)
- **Alignment:** Cell Ranger (10x), STARsolo, Alevin (salmon)
""")

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# NAVIGATION
# ════════════════════════════════════════════════════════════════════════════════

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/10_Public_Data.py",
                 label="← Lesson 10: Public Data & Reproducibility", icon="🌐")
with col_n2:
    st.markdown("*More lessons coming soon.*")

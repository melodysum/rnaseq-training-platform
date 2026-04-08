"""
pages/08_Clustering_Heatmaps.py — Lesson 8: Clustering & Heatmaps
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression, log_cpm
from utils.clustering_utils import (
    cluster_matrix,
    get_leaf_order,
    row_zscore,
    assign_modules,
)
from utils.exploration import sample_correlation_matrix

st.set_page_config(
    page_title="Lesson 8 — Clustering & Heatmaps",
    page_icon="🗺️",
    layout="wide",
)

init_session_data()

st.title("🗺️ Lesson 8 — Clustering & Heatmaps")
st.markdown("""
> **Learning goal:** Understand how to organise RNA-seq expression patterns into
> interpretable clusters of genes and samples, and how heatmaps summarise complex data.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
filter_res = st.session_state.get("filter_results", {})
de_results = st.session_state.get("de_results_lesson7")

MAX_GENES_LIVE = 200

# ── SECTION 1 ─────────────────────────────────────────────────────────────────
st.subheader("🧩 Section 1 — Why clustering matters after DE analysis")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "📋", "DE gives a list",
     "DE tells you *which* genes change — clustering shows *how* they change together across samples."),
    (c2, "🔗", "Clustering reveals patterns",
     "Genes with similar expression profiles may share biological function or regulatory pathway."),
    (c3, "👥", "Sample structure",
     "Sample clustering reveals biological group separation, batch effects, or unexpected outliers."),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:130px;">
<div style="font-size:1.4rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ── SECTION 2 — Sidebar + data loading ───────────────────────────────────────
sources = []
if de_results is not None and de_results["significant"].sum() > 0:
    sources.append("DE significant genes (from Lesson 7)")
sources.append("Top variable genes")
sources.append("Built-in demo selection")

with st.sidebar:
    st.header("🔧 Clustering Controls")
    data_source  = st.selectbox("Gene source", sources)
    n_genes      = st.slider("Number of genes to cluster", 10, 500, 50, 10)

    if n_genes > MAX_GENES_LIVE:
        st.warning(f"⚠️ > {MAX_GENES_LIVE} genes: use **▶ Run Clustering** button.")

    st.divider()
    st.subheader("Transformation")
    use_scaling = st.checkbox("Row Z-score scaling", value=True)

    st.divider()
    st.subheader("Clustering settings")
    dist_metric  = st.selectbox("Distance metric",  ["pearson", "spearman", "euclidean"])
    linkage_meth = st.selectbox("Linkage method",   ["average", "complete", "ward"])
    cluster_genes   = st.checkbox("Cluster genes",   value=True)
    cluster_samples = st.checkbox("Cluster samples", value=True)

    st.divider()
    st.subheader("Gene modules")
    n_modules = st.slider("Number of gene modules (k)", 2, 10, 3)
    run_btn   = st.button("▶ Run Clustering", type="primary")

st.subheader("📥 Section 2 — Input data and gene selection")

counts_use = filter_res.get("counts_filtered", filter_low_expression(counts_raw, 10, 5))
log_expr   = log_cpm(counts_use)

if data_source.startswith("DE significant") and de_results is not None:
    sig_genes  = [g for g in de_results[de_results["significant"]].index if g in log_expr.index]
    gene_pool  = sig_genes
else:
    gene_pool  = log_expr.var(axis=1).nlargest(1000).index.tolist()

selected_genes = gene_pool[:min(n_genes, len(gene_pool))]
expr_subset    = log_expr.loc[selected_genes]

c1, c2, c3, c4 = st.columns(4)
c1.metric("Genes selected", f"{len(selected_genes):,}")
c2.metric("Samples",        f"{expr_subset.shape[1]}")
c3.metric("Source",         data_source.split(" (")[0][:20])
c4.metric("DE results",     "Yes" if de_results is not None else "No")

if len(selected_genes) == 0:
    st.error("No genes available. Try a different source or check that data is loaded.")
    st.stop()

st.divider()

# ── SECTION 3 — Scaling ───────────────────────────────────────────────────────
st.subheader("📏 Section 3 — Transformation and row scaling")

with st.expander("📖 Why row scaling matters", expanded=False):
    col_why, col_ex = st.columns(2)
    with col_why:
        st.markdown("""
**Row Z-score scaling:** `z = (x − mean) / std`

- Each gene is centred at 0
- ±1 = one standard deviation above/below that gene's mean
- All genes visually comparable regardless of absolute expression
- Use a **diverging colour scale** (blue–white–red)

| Goal | Scaling |
|---|---|
| Compare absolute levels | Off |
| Compare up/down patterns | On |
| Most heatmaps in papers | On (Z-score) |
        """)
    with col_ex:
        demo   = pd.DataFrame({"S1":[1000,5,50],"S2":[1200,3,80],"S3":[800,8,30]},
                               index=["HighExpr","LowExpr","MidExpr"])
        demo_z = row_zscore(demo)
        for mat, title, cs, mid in [
            (demo,   "Unscaled", "Viridis", None),
            (demo_z, "Z-scored", "RdBu_r",  0),
        ]:
            kw = dict(z=mat.values, x=list(mat.columns), y=list(mat.index),
                      colorscale=cs, showscale=True)
            if mid is not None:
                kw["zmid"] = mid
            fig = go.Figure(go.Heatmap(**kw))
            fig.update_layout(title=title, height=160, margin=dict(t=30,l=80))
            st.plotly_chart(fig, use_container_width=True)

st.divider()

# ── SECTION 4 — Distance metrics ─────────────────────────────────────────────
st.subheader("📐 Section 4 — Distance metrics and linkage methods")

with st.expander("📖 How metric and linkage choices affect clustering", expanded=False):
    st.markdown("""
| Metric | Measures | Best for |
|---|---|---|
| **Pearson** | Linear profile shape | Co-regulated gene discovery |
| **Spearman** | Rank-based profile shape | Robust, non-linear patterns |
| **Euclidean** | Absolute magnitude difference | Magnitude-sensitive grouping |

| Linkage | Behaviour |
|---|---|
| **Average** | Distance between cluster centres — balanced |
| **Complete** | Maximum pairwise distance — tight clusters |
| **Ward** | Minimises within-cluster variance — often cleanest modules |

Different combinations can give very different results. Try several to check robustness.
    """)

st.divider()

# ── SECTION 5–6 — Clustering + Heatmap ───────────────────────────────────────
st.subheader("🔥 Section 5–6 — Clustering and interactive heatmap")

if use_scaling:
    plot_matrix = row_zscore(expr_subset)
    colorscale  = "RdBu_r"
    zmid        = 0
    color_label = "Z-score"
else:
    plot_matrix = expr_subset.copy()
    colorscale  = "Viridis"
    zmid        = None
    color_label = "log-CPM"

gene_order   = list(plot_matrix.index)
sample_order = list(plot_matrix.columns)
lkg_genes    = None
lkg_samp     = None

should_cluster = run_btn or (n_genes <= MAX_GENES_LIVE)

if should_cluster and (cluster_genes or cluster_samples):
    mat = plot_matrix.values
    if cluster_genes and len(gene_order) > 1:
        try:
            lkg_genes  = cluster_matrix(mat, metric=dist_metric, method=linkage_meth)
            gene_order = [plot_matrix.index[i] for i in get_leaf_order(lkg_genes)]
        except Exception as e:
            st.warning(f"Gene clustering failed: {e}")
            lkg_genes = None
    if cluster_samples and len(sample_order) > 1:
        try:
            lkg_samp    = cluster_matrix(mat.T, metric=dist_metric, method=linkage_meth)
            sample_order = [plot_matrix.columns[i] for i in get_leaf_order(lkg_samp)]
        except Exception as e:
            st.warning(f"Sample clustering failed: {e}")
            lkg_samp = None
elif n_genes > MAX_GENES_LIVE and not run_btn:
    st.info(f"⚡ Clustering paused for >{MAX_GENES_LIVE} genes. Click **▶ Run Clustering** in the sidebar.")

plot_ordered = plot_matrix.loc[gene_order, sample_order]

hover_text = []
for gene in gene_order:
    row = []
    for sample in sample_order:
        raw_val = expr_subset.loc[gene, sample]
        z_val   = plot_matrix.loc[gene, sample]
        if use_scaling:
            row.append(f"Gene: {gene}<br>Sample: {sample}<br>log-CPM: {raw_val:.2f}<br>Z-score: {z_val:.2f}")
        else:
            row.append(f"Gene: {gene}<br>Sample: {sample}<br>log-CPM: {raw_val:.2f}")
    hover_text.append(row)

hm_kw = dict(
    z=plot_ordered.values, x=sample_order, y=gene_order,
    colorscale=colorscale, text=hover_text, hoverinfo="text",
    colorbar=dict(title=color_label, len=0.8),
)
if zmid is not None:
    hm_kw["zmid"] = zmid

fig_hm = go.Figure(go.Heatmap(**hm_kw))
fig_hm.update_layout(
    height=max(400, min(900, len(gene_order) * 14 + 100)),
    margin=dict(l=130, r=60, t=50, b=110),
    xaxis=dict(tickangle=45, tickfont=dict(size=9 if len(sample_order) > 20 else 11)),
    yaxis=dict(tickfont=dict(size=8 if len(gene_order) > 50 else 10),
               showticklabels=(len(gene_order) <= 60)),
    title=f"{'Z-scored' if use_scaling else 'Unscaled'} heatmap — "
          f"{len(gene_order)} genes × {len(sample_order)} samples",
)
st.plotly_chart(fig_hm, use_container_width=True)

if len(gene_order) > 60:
    st.caption("Gene labels hidden for readability. Hover over cells to see gene IDs.")

meta_cols = [c for c in ["groupA", "batch"] if c in metadata.columns]
if meta_cols:
    with st.expander("Sample annotation"):
        ann_samples = [s for s in sample_order if s in metadata.index]
        st.dataframe(metadata.loc[ann_samples, meta_cols].T, use_container_width=True)

st.divider()

# ── SECTION 7 ─────────────────────────────────────────────────────────────────
st.subheader("🔄 Section 7 — Sample clustering vs gene clustering")

st.markdown("""
| Mode | What it reveals |
|---|---|
| **Gene clustering ON** | Co-expressed gene groups — potential co-regulation or shared function |
| **Sample clustering ON** | Cohort structure — biological groups, batch effects, outliers |
| **Both ON** | Block structure: which gene groups are up/down in which sample groups |
| **Both OFF** | Original order — useful as a reference comparison |

Toggle these in the sidebar to explore how the heatmap changes.
""")

st.divider()

# ── SECTION 8 — Gene modules ──────────────────────────────────────────────────
st.subheader("🧬 Section 8 — Identifying gene modules")

if lkg_genes is not None and len(gene_order) > 1:
    module_assignments = assign_modules(lkg_genes, n_modules, gene_order)
    mod_summary = module_assignments.value_counts().sort_index().reset_index()
    mod_summary.columns = ["Module", "Gene count"]

    col_a, col_b = st.columns([2, 3])
    with col_a:
        fig_mod = px.bar(
            mod_summary, x="Module", y="Gene count",
            color="Module", title=f"{n_modules} gene modules",
            height=280,
        )
        fig_mod.update_layout(margin=dict(t=40), showlegend=False)
        st.plotly_chart(fig_mod, use_container_width=True)

    with col_b:
        sel_mod   = st.selectbox(
            "Inspect module",
            sorted(module_assignments.unique()),
            format_func=lambda x: f"Module {x} ({(module_assignments==x).sum()} genes)",
        )
        mod_genes = module_assignments[module_assignments == sel_mod].index.tolist()
        st.markdown(f"**Module {sel_mod} — {len(mod_genes)} genes:**")
        st.dataframe(pd.DataFrame({"Gene": mod_genes}),
                     use_container_width=True, height=200, hide_index=True)
        st.download_button(
            f"⬇️ Module {sel_mod} gene list (.txt)",
            data="\n".join(mod_genes),
            file_name=f"module_{sel_mod}_genes.txt",
            mime="text/plain",
        )

    st.download_button(
        "⬇️ All module assignments (CSV)",
        data=module_assignments.reset_index().rename(columns={"index":"gene"}).to_csv(index=False),
        file_name="gene_module_assignments.csv",
        mime="text/csv",
    )
    st.info("Module gene lists can be used as input for pathway enrichment analysis in Lesson 9.")
else:
    st.info("Enable **Cluster genes** and click **▶ Run Clustering** to identify modules.")

st.divider()

# ── SECTION 9 — Sample correlation ────────────────────────────────────────────
st.subheader("🗺️ Section 9 — Sample-to-sample correlation")

corr_mat = sample_correlation_matrix(log_expr)
corr_ord = corr_mat.loc[sample_order, sample_order]
fig_corr = px.imshow(
    corr_ord, color_continuous_scale="RdBu",
    zmin=0.8, zmax=1.0,
    title="Sample Pearson correlation (log-CPM, all filtered genes)",
    height=460, labels=dict(color="r"),
)
fig_corr.update_layout(margin=dict(t=40))
st.plotly_chart(fig_corr, use_container_width=True)
st.caption("If samples cluster by batch rather than by biology, revisit Lesson 6 — Batch Correction.")

st.divider()

# ── SECTION 10 — Common mistakes ──────────────────────────────────────────────
st.subheader("⚠️ Section 10 — Common interpretation mistakes")

st.error("**Too many genes → unreadable heatmap.** Stick to 50–200 meaningful genes.")
st.warning("**Forgetting that scaling changes interpretation.** Scaled and unscaled heatmaps are not comparable.")
st.warning("**Overinterpreting colour blocks.** Red = relatively high for that gene, not necessarily biologically important.")
st.info("**Assuming one result is the truth.** Different metrics, linkage, and k can all give different patterns.")
st.info("**Clustering random genes.** Always start from DE-significant or highly variable genes.")

st.divider()

# ── SECTION 11 — Takeaways ────────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

cols = st.columns(5)
for col, (icon, title, body) in zip(cols, [
    ("🗺️", "Patterns after DE", "Clustering organises gene lists into interpretable expression programmes."),
    ("📏", "Z-score scaling", "Centres each gene to highlight relative up/down across samples."),
    ("📐", "Metric matters", "Pearson/Spearman focus on shape; Euclidean on magnitude."),
    ("🧬", "Modules → enrichment", "Gene modules can feed directly into pathway analysis in Lesson 9."),
    ("🔍", "Exploratory only", "Clustering generates hypotheses — it is not a statistical test."),
]):
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;font-size:0.9rem;margin:0.3rem 0">{title}</div>
<div style="font-size:0.82rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/07_Differential_Expression.py",
                 label="← Lesson 7: Differential Expression", icon="🧪")
with col_n2:
    st.page_link("pages/09_Functional_Enrichment.py",
                 label="Lesson 9: Functional Enrichment →", icon="🧩")

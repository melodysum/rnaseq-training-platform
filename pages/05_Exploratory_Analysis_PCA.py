"""
pages/05_Exploratory_Analysis_PCA.py — Lesson 4: Exploratory Analysis & PCA
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression, log_cpm
from utils.batch_effects import top_variable_genes
from utils.pca_utils import run_pca, pca_plot_df
from utils.exploration import (
    sample_correlation_matrix,
    sample_distance_matrix,
    outlier_scores,
)
from utils.batch_effects import top_variable_genes as tvg

st.set_page_config(
    page_title="Lesson 5 — Exploratory Analysis",
    page_icon="📊",
    layout="wide",
)

init_session_data()

st.title("📊 Lesson 5 — Exploratory Analysis & PCA")
st.markdown("""
> **Learning goal:** Learn how to inspect RNA-seq sample structure, detect outliers
> and batch effects, and correctly interpret PCA before running differential expression.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]

# ── Decide which count matrix to use ─────────────────────────────────────────
filter_res = st.session_state.get("filter_results", {})

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔧 Analysis Controls")

    use_filtered = st.radio(
        "Gene set",
        ["Use filtered genes (from Lesson 1)", "Use all genes"],
        index=0 if filter_res else 1,
    )
    if use_filtered == "Use filtered genes (from Lesson 1)" and filter_res:
        counts_use = filter_res["counts_filtered"]
        st.caption(f"Using {len(counts_use):,} filtered genes.")
    else:
        counts_use = filter_low_expression(counts_raw, 10, 5)
        st.caption(f"Using {len(counts_use):,} genes (default filter).")

    st.divider()
    top_n_genes = st.slider("Top variable genes for PCA", 100,
                             min(3000, len(counts_use)), 500, 100)
    pc_choices  = [f"PC{i}" for i in range(1, 6)]
    pc_x = st.selectbox("PCA X axis", pc_choices, index=0)
    pc_y = st.selectbox("PCA Y axis", pc_choices, index=1)

    meta_cols = [c for c in metadata.columns
                 if metadata[c].nunique() <= 20 or c == "age"]
    color_by = st.selectbox("Colour samples by", meta_cols)
    shape_by = st.selectbox("Shape samples by (optional)",
                             ["None"] + [c for c in meta_cols if c != color_by])
    show_labels = st.checkbox("Show sample labels", value=False)

# ── Preprocess ────────────────────────────────────────────────────────────────
log_expr = log_cpm(counts_use)
scores, explained, loadings = run_pca(log_expr, n_components=5,
                                       top_var_genes=top_n_genes)
pca_df = pca_plot_df(scores, metadata, pc_x, pc_y)

# ── SECTION 1 — Why EDA matters ───────────────────────────────────────────────
st.subheader("🔭 Section 1 — Why exploratory analysis matters")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "🕵️", "Detect outliers",
     "Samples that look very different from all others may have quality issues, contamination, or unexpected biology."),
    (c2, "🔄", "Reveal batch effects",
     "If samples separate by sequencing run or library prep date rather than biology, you need to account for this."),
    (c3, "✅", "Confirm study design",
     "PCA confirms that your groups separate as expected — giving you confidence before running DE analysis."),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:150px;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ── SECTION 2 — Dataset overview ─────────────────────────────────────────────
st.subheader("📋 Section 2 — Dataset overview")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown(f"""
| | |
|---|---|
| **Genes used** | {len(counts_use):,} |
| **Samples** | {len(counts_use.columns)} |
| **Metadata columns** | {', '.join(metadata.columns)} |
| **Groups** | {', '.join(metadata['groupA'].unique())} |
| **Filtering applied** | {'Yes (from Lesson 1)' if filter_res and use_filtered.startswith('Use filtered') else 'Default (min 10 counts, 5 samples)'} |
""")
with col_r:
    st.info("""
**Before running DE analysis, always check:**
- Do samples of the same group cluster together?
- Is there any unexpected separation by a technical variable?
- Are there any obvious outlier samples?
    """)

st.divider()

# ── SECTION 3 — Transformation ────────────────────────────────────────────────
st.subheader("🔢 Section 3 — Why we transform counts before PCA")

col_why, col_effect = st.columns(2)
with col_why:
    st.markdown("""
**Raw counts are not suitable for PCA because:**
- Library sizes differ between samples — a gene with 1000 counts in a
  5M-read library expresses differently than in a 20M-read library
- Highly expressed genes dominate PCA if counts are used directly
- The variance of count data scales with the mean (heteroscedasticity)

**log₂(CPM + 1) solves this by:**
- Normalising for library size (CPM step)
- Compressing the dynamic range (log step)
- Stabilising variance across expression levels

**Focusing on top variable genes:**
- Keeps the most informative genes
- Reduces noise from unexpressed/constitutive genes
- Makes PCA faster and cleaner
    """)
with col_effect:
    # Show variance vs mean for raw vs logCPM
    gene_means_raw = counts_use.mean(axis=1)
    gene_vars_raw  = counts_use.var(axis=1)
    gene_means_log = log_expr.mean(axis=1)
    gene_vars_log  = log_expr.var(axis=1)

    sample_idx = np.random.default_rng(42).choice(
        len(counts_use), size=min(1500, len(counts_use)), replace=False)

    fig_mv = go.Figure()
    fig_mv.add_trace(go.Scatter(
        x=np.log2(gene_means_raw.iloc[sample_idx] + 1),
        y=np.log2(gene_vars_raw.iloc[sample_idx] + 1),
        mode="markers", name="Raw counts",
        marker=dict(color="#94a3b8", size=3, opacity=0.5),
    ))
    fig_mv.add_trace(go.Scatter(
        x=gene_means_log.iloc[sample_idx],
        y=gene_vars_log.iloc[sample_idx],
        mode="markers", name="log-CPM",
        marker=dict(color="#3b82f6", size=3, opacity=0.5),
    ))
    fig_mv.update_layout(
        xaxis_title="Mean expression (log scale)",
        yaxis_title="Variance (log scale)",
        height=300, margin=dict(t=10),
        legend=dict(x=0.6, y=0.95),
        title="Mean-variance relationship",
    )
    st.plotly_chart(fig_mv, use_container_width=True)
    st.caption("log-CPM (blue) has a flatter mean-variance trend — better for PCA.")

st.divider()

# ── SECTION 4 — PCA ───────────────────────────────────────────────────────────
st.subheader("🔵 Section 4 — PCA analysis")

col_pca, col_scree = st.columns([3, 2])

with col_pca:
    shape_col = None if shape_by == "None" else shape_by
    fig_pca = px.scatter(
        pca_df, x=pc_x, y=pc_y,
        color=color_by,
        symbol=shape_col,
        text="sample" if show_labels else None,
        labels={
            pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
            pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)",
        },
        title=f"PCA — coloured by {color_by}",
        height=430,
    )
    fig_pca.update_traces(marker=dict(size=10), textposition="top center")
    fig_pca.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_pca, use_container_width=True)

with col_scree:
    exp_df = pd.DataFrame({
        "PC":   [f"PC{i+1}" for i in range(len(explained))],
        "Var %": explained * 100,
        "Cum %": np.cumsum(explained) * 100,
    })
    fig_scree = go.Figure()
    fig_scree.add_trace(go.Bar(
        x=exp_df["PC"], y=exp_df["Var %"],
        name="Per PC", marker_color="#3b82f6",
    ))
    fig_scree.add_trace(go.Scatter(
        x=exp_df["PC"], y=exp_df["Cum %"],
        name="Cumulative", line=dict(color="#ef4444"),
        yaxis="y2", mode="lines+markers",
    ))
    fig_scree.update_layout(
        title="Scree plot",
        yaxis=dict(title="Variance explained (%)"),
        yaxis2=dict(title="Cumulative (%)", overlaying="y",
                    side="right", range=[0, 100]),
        height=430, margin=dict(t=40),
        legend=dict(x=0.5, y=0.95),
    )
    st.plotly_chart(fig_scree, use_container_width=True)

    st.dataframe(
        exp_df.style.format({"Var %": "{:.2f}", "Cum %": "{:.2f}"}),
        use_container_width=True, hide_index=True, height=180,
    )

with st.expander("💡 How to interpret this PCA"):
    st.markdown(f"""
- **{pc_x}** explains **{explained[int(pc_x[-1])-1]*100:.1f}%** of variance —
  the largest single source of expression variation in this dataset.
- **{pc_y}** explains **{explained[int(pc_y[-1])-1]*100:.1f}%** — the second largest.
- Samples close together have **more similar expression profiles**.
- PCA does not automatically label axes as "biology" or "batch" —
  you need to overlay metadata to interpret what the separation means.
- Switch the colour variable in the sidebar to see how different metadata
  variables align with the PC structure.
    """)

st.divider()

# ── SECTION 5 — Outlier detection ────────────────────────────────────────────
st.subheader("🚨 Section 5 — Outlier detection")

outlier_z = outlier_scores(scores, pcs=[pc_x, pc_y])
pca_df["outlier_score"] = pca_df["sample"].map(outlier_z)
threshold = 2.0

col_out1, col_out2 = st.columns(2)

with col_out1:
    pca_df["outlier"] = pca_df["outlier_score"] > threshold
    fig_out = px.scatter(
        pca_df, x=pc_x, y=pc_y,
        color="outlier",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        hover_name="sample",
        hover_data={col: True for col in metadata.columns
                    if col in pca_df.columns},
        labels={
            pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
            pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)",
        },
        title="Potential outliers (z-score > 2 from centroid)",
        height=370,
    )
    fig_out.update_traces(marker=dict(size=9))
    fig_out.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_out, use_container_width=True)

with col_out2:
    outlier_df = pca_df[["sample", "outlier_score", "outlier"]].copy()
    if color_by in pca_df.columns:
        outlier_df[color_by] = pca_df[color_by].values
    outlier_df = outlier_df.sort_values("outlier_score", ascending=False)

    st.markdown("**Outlier scores (all samples)**")
    st.dataframe(
        outlier_df.style
        .format({"outlier_score": "{:.2f}"})
        .apply(lambda col: [
            "background-color: #fee2e2" if v else ""
            for v in col
        ], subset=["outlier"]),
        use_container_width=True, hide_index=True, height=320,
    )

    n_outliers = outlier_df["outlier"].sum()
    if n_outliers > 0:
        st.warning(
            f"**{n_outliers} potential outlier(s) flagged.** "
            "Inspect their metadata and raw expression before deciding to exclude them."
        )
    else:
        st.success("No obvious outliers detected at z-score > 2.")

st.caption("""
⚠️ Outlier detection here is based on PCA distance from the centroid.
This is an exploratory flag, not a definitive removal criterion.
Always check the biology and QC metrics before excluding any sample.
""")

st.divider()

# ── SECTION 6 — Sample distance heatmap ──────────────────────────────────────
st.subheader("🗺️ Section 6 — Sample-to-sample distances")

dist_mat = sample_distance_matrix(log_expr)
corr_mat = sample_correlation_matrix(log_expr)

tab_dist, tab_corr = st.tabs(["Distance heatmap", "Correlation heatmap"])

with tab_dist:
    fig_dist = px.imshow(
        dist_mat,
        color_continuous_scale="Blues_r",
        title="Euclidean distance between samples (log-CPM)",
        height=500,
        labels=dict(color="Distance"),
    )
    fig_dist.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_dist, use_container_width=True)
    st.caption(
        "Darker = more similar. Samples from the same group should generally "
        "be darker (more similar) to each other than to samples from other groups."
    )

with tab_corr:
    fig_corr = px.imshow(
        corr_mat,
        color_continuous_scale="RdBu",
        zmin=0.8, zmax=1.0,
        title="Pearson correlation between samples (log-CPM)",
        height=500,
        labels=dict(color="r"),
    )
    fig_corr.update_layout(margin=dict(t=40))
    st.plotly_chart(fig_corr, use_container_width=True)
    st.caption(
        "Values close to 1.0 indicate highly similar samples. "
        "Low-correlation samples (unusual rows/columns) may be outliers."
    )

st.divider()

# ── SECTION 7 — Metadata interpretation ──────────────────────────────────────
st.subheader("🏷️ Section 7 — Interpreting metadata on PCA")

st.markdown("""
Switch through the metadata columns below to understand what drives sample separation.
""")

meta_numeric = [c for c in metadata.columns
                if pd.api.types.is_numeric_dtype(metadata[c])]
meta_categ   = [c for c in metadata.columns
                if not pd.api.types.is_numeric_dtype(metadata[c])]

n_cols   = min(len(meta_cols), 4)
tab_list = st.tabs(meta_cols[:n_cols])

for tab, col_name in zip(tab_list, meta_cols[:n_cols]):
    with tab:
        pca_tab = pca_plot_df(scores, metadata, pc_x, pc_y)
        color_seq = px.colors.qualitative.Set2

        if col_name in meta_numeric:
            fig_t = px.scatter(
                pca_tab, x=pc_x, y=pc_y, color=col_name,
                color_continuous_scale="Viridis",
                hover_name="sample", height=360,
                labels={pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
                        pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)"},
                title=f"Coloured by {col_name}",
            )
        else:
            fig_t = px.scatter(
                pca_tab, x=pc_x, y=pc_y, color=col_name,
                color_discrete_sequence=color_seq,
                hover_name="sample", height=360,
                labels={pc_x: f"{pc_x} ({explained[int(pc_x[-1])-1]*100:.1f}%)",
                        pc_y: f"{pc_y} ({explained[int(pc_y[-1])-1]*100:.1f}%)"},
                title=f"Coloured by {col_name}",
            )
        fig_t.update_traces(marker=dict(size=10))
        fig_t.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_t, use_container_width=True)

        guidance = {
            "groupA": "✅ If groups separate here, that supports a real biological difference.",
            "batch":  "⚠️ If batches separate, consider batch correction before DE analysis.",
            "donor":  "ℹ️ Donor effects are common in paired designs — use a paired test.",
            "sex":    "ℹ️ Sex differences can be a real biological signal — check if relevant.",
            "age":    "ℹ️ Age gradients can reflect biology — worth noting in your model.",
        }
        msg = guidance.get(col_name, f"ℹ️ Consider whether {col_name} should be modelled.")
        st.caption(msg)

st.divider()

# ── SECTION 8 — Common mistakes ───────────────────────────────────────────────
st.subheader("⚠️ Section 8 — Common PCA interpretation mistakes")

st.error("**Assuming separation always means biology.** Batch effects, library size differences, or sample swaps can all cause PCA separation.")
st.warning("**Ignoring metadata.** PCA axes have no inherent biological meaning until you overlay metadata.")
st.warning("**Overinterpreting weak PCs.** PC3 explaining 3% of variance may just be noise.")
st.info("**Running PCA on raw counts.** Always transform first — raw counts give misleading PCA due to library size effects.")
st.info("**Treating PCA as a statistical test.** PCA is exploratory — it generates hypotheses, not conclusions.")

st.divider()

# ── SECTION 9 — Takeaways ─────────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

cols = st.columns(4)
msgs = [
    ("🔭", "EDA first", "Always explore your data before running DE analysis."),
    ("🔢", "Transform counts", "log-CPM stabilises variance for PCA."),
    ("🏷️", "Use metadata", "PCA interpretation requires metadata context."),
    ("🚨", "Check outliers", "Identify unusual samples before DE — don't ignore them."),
]
for col, (icon, title, body) in zip(cols, msgs):
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;">
<div style="font-size:1.5rem">{icon}</div>
<div style="font-weight:600;font-size:0.9rem;margin:0.3rem 0">{title}</div>
<div style="font-size:0.82rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — VARIANCE DECOMPOSITION (PVCA-like)
# ══════════════════════════════════════════════════════════════════════════════

from utils.pca_utils import variance_decomposition

st.subheader("📊 Section 9 — Variance Decomposition")
st.markdown("""
PCA plots show *whether* samples cluster by batch or condition.
Variance decomposition quantifies *how much* of total expression variance
each factor explains — giving a number, not just a visual impression.

**Method (PVCA-like, Bushel et al. 2009):**
1. Select PCs that together explain ≥ 80% of variance.
2. For each PC, regress scores against experimental factors (OLS).
3. Compute partial R² for each factor within each PC.
4. Weight by each PC's share of variance and aggregate.

**Example interpretation:**
- Batch explains 35%, condition explains 40% → batch correction is warranted.
- Batch explains 8%, condition explains 55% → batch is minor, proceed.
""")

with st.expander("⚙️ Configure decomposition", expanded=True):
    vd_factors = []
    meta_cols_vd = [c for c in metadata.columns if metadata[c].nunique() < 20]
    if meta_cols_vd:
        vd_factors = st.multiselect(
            "Factors to decompose",
            options=meta_cols_vd,
            default=meta_cols_vd[:3] if len(meta_cols_vd) >= 3 else meta_cols_vd,
            help="Select metadata columns. Categorical and continuous both work.",
        )
    vd_cumvar = st.slider(
        "Include PCs up to cumulative variance (%)",
        50, 95, 80, 5,
        key="vd_cumvar",
    )

if vd_factors:
    with st.spinner("Running variance decomposition…"):
        try:
            vd_result = variance_decomposition(
                scores, explained,
                metadata,
                factors=vd_factors,
                cumvar_threshold=vd_cumvar / 100,
            )

            col_vd1, col_vd2 = st.columns([1, 1])
            with col_vd1:
                st.dataframe(
                    vd_result.style.format({"pct_variance_explained": "{:.1f}%"}),
                    use_container_width=True,
                    hide_index=True,
                )
            with col_vd2:
                import plotly.express as px
                fig_vd = px.bar(
                    vd_result,
                    x="factor",
                    y="pct_variance_explained",
                    color="factor",
                    labels={"pct_variance_explained": "% variance explained",
                            "factor": "Factor"},
                    color_discrete_sequence=["#2563eb", "#16a34a", "#dc2626",
                                             "#d97706", "#7c3aed", "#94a3b8"],
                )
                fig_vd.update_layout(
                    height=350,
                    showlegend=False,
                    yaxis=dict(range=[0, 100]),
                )
                st.plotly_chart(fig_vd, use_container_width=True)

            top = vd_result[vd_result["factor"] != "Residual / other"].iloc[0]
            st.info(
                f"**Dominant factor:** {top['factor']} explains "
                f"{top['pct_variance_explained']:.1f}% of total variance. "
                "Values above ~30% for batch indicate batch correction is needed "
                "before differential expression analysis."
            )
        except Exception as e:
            st.error(f"Variance decomposition failed: {e}")
else:
    st.info("Select at least one factor above to run decomposition.")

st.divider()
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/06_Batch_Correction.py", label="← Lesson 6: Batch Correction", icon="🔄")
with col_n2:
    st.page_link("pages/07_Differential_Expression.py", label="Lesson 7: Differential Expression →", icon="🧪")

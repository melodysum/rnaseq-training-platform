"""
pages/06_Batch_Correction.py — Lesson 3: Batch Effects & Correction
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression
from utils.batch_effects import (
    to_log_cpm,
    simple_batch_correction,
    variance_explained_by,
    top_variable_genes,
)
from utils.pca_utils import run_pca, pca_plot_df
from utils.simulation import simulate_batch_data

st.set_page_config(
    page_title="Lesson 6 — Batch Correction",
    page_icon="🔄",
    layout="wide",
)

init_session_data()

# ── Header ────────────────────────────────────────────────────────────────────
st.title("🔄 Lesson 6 — Batch Effects & Correction")
st.markdown("""
> **Learning goal:** Understand how technical batch effects arise, how to detect them
> with PCA, how to correct for them, and — crucially — when correction can go wrong.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
has_batch  = "batch" in metadata.columns

# ── Preprocessing shared across sections ──────────────────────────────────────
filter_res = st.session_state.get("filter_results", {})
if filter_res:
    counts_use = filter_res["counts_filtered"]
else:
    counts_use = filter_low_expression(counts_raw, 10, 5)

log_expr = to_log_cpm(counts_use)

# ── SECTION 1 — What is a batch effect? ──────────────────────────────────────
st.subheader("🧩 Section 1 — What is a batch effect?")

col_def, col_diag = st.columns([3, 2])
with col_def:
    st.markdown("""
**Batch effects** are unwanted technical differences between groups of samples
that are unrelated to the biology you care about.

They arise from things like:
- Different sequencing runs or machines
- Different operators or reagent lots
- Different library preparation dates
- Different labs or institutions

The danger is that batch effects can make samples **separate by technical variables**
rather than by biology — confusing your downstream analysis.
    """)

with col_diag:
    # Simple schematic using a plotly figure
    fig_schematic = go.Figure()
    labels = ["Biological signal", "Batch effect", "What you observe"]
    colors = ["#22c55e", "#f87171", "#3b82f6"]
    x_pos  = [0.15, 0.50, 0.85]
    for x, label, color in zip(x_pos, labels, colors):
        fig_schematic.add_shape(
            type="rect", x0=x-0.12, x1=x+0.12, y0=0.3, y1=0.7,
            fillcolor=color, opacity=0.85, line_width=0,
        )
        fig_schematic.add_annotation(
            x=x, y=0.5, text=label, showarrow=False,
            font=dict(color="white", size=11, family="Inter"),
            align="center",
        )
    # Arrows
    fig_schematic.add_annotation(
        x=0.48, y=0.5, ax=0.27, ay=0.5,
        arrowhead=2, arrowcolor="#64748b", arrowwidth=2,
        text="+", font=dict(size=16, color="#64748b"), showarrow=True,
    )
    fig_schematic.add_annotation(
        x=0.73, y=0.5, ax=0.62, ay=0.5,
        arrowhead=2, arrowcolor="#64748b", arrowwidth=2,
        text="=", font=dict(size=16, color="#64748b"), showarrow=True,
    )
    fig_schematic.update_layout(
        height=200, margin=dict(l=0, r=0, t=10, b=10),
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[0, 1]),
        plot_bgcolor="white",
    )
    st.plotly_chart(fig_schematic, use_container_width=True)
    st.caption("Observed expression = true biology + technical batch shift + noise")

st.divider()

# ── SECTION 2 — Dataset overview ─────────────────────────────────────────────
st.subheader("📋 Section 2 — Dataset overview")

col_meta, col_warn = st.columns([2, 2])
with col_meta:
    n_genes   = len(counts_use)
    n_samples = len(counts_use.columns)
    meta_cols = list(metadata.columns)
    st.markdown(f"""
| | |
|---|---|
| **Genes (after filtering)** | {n_genes:,} |
| **Samples** | {n_samples} |
| **Metadata columns** | {', '.join(meta_cols)} |
| **Batch column detected** | {'✅ Yes' if has_batch else '❌ No'} |
""")

with col_warn:
    if has_batch:
        batch_counts = metadata["batch"].value_counts()
        st.success(
            f"✅ Batch column found with **{len(batch_counts)} batches**: "
            + ", ".join(f"{b} ({n} samples)" for b, n in batch_counts.items())
        )
    else:
        st.warning(
            "⚠️ No `batch` column found in your metadata. "
            "The correction demonstration will use simulated data. "
            "Upload metadata with a `batch` column to use your own data."
        )

st.divider()

# ── SECTION 3 — PCA: detecting batch effects ──────────────────────────────────
st.subheader("🔭 Section 3 — Detecting batch effects with PCA")

st.markdown("""
PCA is the first tool to reach for when checking for batch effects.
If samples cluster more strongly by batch than by biology, that's a warning sign.
""")

with st.sidebar:
    st.header("🔧 PCA Controls")
    top_genes = st.slider("Top variable genes for PCA", 100, min(2000, n_genes), 500, 100)
    pc_x = st.selectbox("X axis", ["PC1", "PC2", "PC3"], index=0)
    pc_y = st.selectbox("Y axis", ["PC1", "PC2", "PC3"], index=1)
    show_labels = st.checkbox("Show sample labels", value=False)
    color_by = st.selectbox(
        "Colour by",
        [c for c in ["groupA", "batch", "sex", "donor"] if c in metadata.columns],
    )

scores, explained, _ = run_pca(log_expr, n_components=5, top_var_genes=top_genes)
pca_df = pca_plot_df(scores, metadata, pc_x, pc_y)

col_bio, col_batch = st.columns(2)

def make_pca_fig(df, color_col, title, explained, pc_x, pc_y, show_labels):
    x_var = explained[int(pc_x[-1]) - 1] * 100
    y_var = explained[int(pc_y[-1]) - 1] * 100
    fig = px.scatter(
        df, x=pc_x, y=pc_y,
        color=color_col if color_col in df.columns else None,
        text="sample" if show_labels else None,
        title=title,
        labels={
            pc_x: f"{pc_x} ({x_var:.1f}%)",
            pc_y: f"{pc_y} ({y_var:.1f}%)",
        },
        height=380,
    )
    fig.update_traces(marker=dict(size=9), textposition="top center")
    fig.update_layout(margin=dict(t=40))
    return fig

with col_bio:
    fig_bio = make_pca_fig(pca_df, "groupA", "PCA coloured by biological group",
                           explained, pc_x, pc_y, show_labels)
    st.plotly_chart(fig_bio, use_container_width=True)

with col_batch:
    color_col = "batch" if has_batch else "groupA"
    title_str = "PCA coloured by batch" if has_batch else "PCA (no batch column — coloured by group)"
    fig_bat = make_pca_fig(pca_df, color_col, title_str,
                           explained, pc_x, pc_y, show_labels)
    st.plotly_chart(fig_bat, use_container_width=True)

# Variance explained bar chart
exp_df = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(len(explained))],
    "Variance explained (%)": explained * 100,
})
fig_var = px.bar(
    exp_df, x="PC", y="Variance explained (%)",
    title="Variance explained per PC",
    color_discrete_sequence=["#3b82f6"],
    height=260,
)
fig_var.update_layout(margin=dict(t=40, b=10))
st.plotly_chart(fig_var, use_container_width=True)

if has_batch:
    r2_bio   = variance_explained_by(scores, metadata, "groupA")
    r2_batch = variance_explained_by(scores, metadata, "batch")
    st.info(
        f"PC1+PC2 variance explained: "
        f"**biological group ≈ {r2_bio:.1%}** | "
        f"**batch ≈ {r2_batch:.1%}**. "
        + ("Batch explains more — correction may help." if r2_batch > r2_bio
           else "Biology dominates — batch effect appears modest.")
    )

st.divider()

# ── SECTION 4 — Simulation ────────────────────────────────────────────────────
st.subheader("🧪 Section 4 — Simulating batch effects")

st.markdown("""
Use the controls below to simulate expression data with varying batch effect strengths.
Compare a **balanced design** (each group appears in all batches) vs a **confounded design**
(each group maps to one batch).
""")

sim_col, sim_out = st.columns([1, 2])

with sim_col:
    sim_n_genes   = st.slider("Genes", 100, 1000, 300, 100, key="sim_genes")
    sim_n_samp    = st.slider("Samples per group", 3, 10, 5, 1, key="sim_samp")
    sim_bio       = st.slider("Biological effect strength", 0.0, 5.0, 2.0, 0.5, key="sim_bio")
    sim_batch     = st.slider("Batch effect strength",      0.0, 5.0, 2.0, 0.5, key="sim_bat")
    sim_confound  = st.checkbox("Confounded design (dangerous!)", value=False, key="sim_conf")

sim_expr, sim_meta = simulate_batch_data(
    n_genes=sim_n_genes,
    n_samples_per_group=sim_n_samp,
    bio_effect=sim_bio,
    batch_effect=sim_batch,
    confounded=sim_confound,
)

sim_scores, sim_exp, _ = run_pca(sim_expr, n_components=3, scale=True)
sim_pca_df = pca_plot_df(sim_scores, sim_meta, "PC1", "PC2")

with sim_out:
    tab_bio, tab_bat = st.tabs(["Colour by biology", "Colour by batch"])
    with tab_bio:
        fig_s1 = px.scatter(
            sim_pca_df, x="PC1", y="PC2", color="groupA",
            title="Simulated PCA — biological groups",
            labels={"PC1": f"PC1 ({sim_exp[0]*100:.1f}%)",
                    "PC2": f"PC2 ({sim_exp[1]*100:.1f}%)"},
            height=350,
        )
        fig_s1.update_traces(marker=dict(size=10))
        st.plotly_chart(fig_s1, use_container_width=True)
    with tab_bat:
        fig_s2 = px.scatter(
            sim_pca_df, x="PC1", y="PC2", color="batch",
            title="Simulated PCA — batch",
            labels={"PC1": f"PC1 ({sim_exp[0]*100:.1f}%)",
                    "PC2": f"PC2 ({sim_exp[1]*100:.1f}%)"},
            color_discrete_sequence=["#f97316", "#8b5cf6"],
            height=350,
        )
        fig_s2.update_traces(marker=dict(size=10))
        st.plotly_chart(fig_s2, use_container_width=True)

    if sim_confound:
        st.error(
            "⚠️ **Confounded design**: batch and biology align perfectly. "
            "Any correction here will also remove real biological signal. "
            "This is why study design matters more than any computational fix."
        )
    else:
        st.success(
            "✅ **Balanced design**: batches are spread across biological groups. "
            "Correction can separate technical from biological variation."
        )

st.divider()

# ── SECTION 5 — Before vs after correction ───────────────────────────────────
st.subheader("⚖️ Section 5 — Before vs after batch correction")

if has_batch:
    corrected_expr = simple_batch_correction(log_expr, metadata, "batch", "groupA")
    data_label = "real uploaded data"
else:
    st.info(
        "No batch column in your metadata — using simulated data for the correction demo."
    )
    sim_demo, meta_demo = simulate_batch_data(
        n_genes=500, n_samples_per_group=5,
        bio_effect=2.0, batch_effect=3.0, confounded=False,
    )
    log_expr_demo   = sim_demo
    corrected_expr  = simple_batch_correction(sim_demo, meta_demo, "batch", "groupA")
    metadata        = meta_demo
    log_expr        = log_expr_demo
    has_batch       = True
    data_label      = "simulated demo data"

scores_before, exp_before, _ = run_pca(
    top_variable_genes(log_expr, 500), n_components=3)
scores_after,  exp_after,  _ = run_pca(
    top_variable_genes(corrected_expr, 500), n_components=3)

pca_before = pca_plot_df(scores_before, metadata, "PC1", "PC2")
pca_after  = pca_plot_df(scores_after,  metadata, "PC1", "PC2")

col_before, col_after = st.columns(2)

with col_before:
    st.markdown("**Before correction**")
    for color_var, tab_label in [("groupA", "Group"), ("batch", "Batch")]:
        if color_var not in pca_before.columns:
            continue
        fig_b = px.scatter(
            pca_before, x="PC1", y="PC2", color=color_var,
            labels={"PC1": f"PC1 ({exp_before[0]*100:.1f}%)",
                    "PC2": f"PC2 ({exp_before[1]*100:.1f}%)"},
            height=320, title=f"Coloured by {color_var}",
        )
        fig_b.update_traces(marker=dict(size=9))
        fig_b.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_b, use_container_width=True)

with col_after:
    st.markdown("**After correction**")
    for color_var in ["groupA", "batch"]:
        if color_var not in pca_after.columns:
            continue
        fig_a = px.scatter(
            pca_after, x="PC1", y="PC2", color=color_var,
            labels={"PC1": f"PC1 ({exp_after[0]*100:.1f}%)",
                    "PC2": f"PC2 ({exp_after[1]*100:.1f}%)"},
            height=320, title=f"Coloured by {color_var}",
        )
        fig_a.update_traces(marker=dict(size=9))
        fig_a.update_layout(margin=dict(t=40))
        st.plotly_chart(fig_a, use_container_width=True)

st.caption(
    f"Using {data_label}. "
    "Correction should reduce batch clustering while preserving biological separation."
)

st.divider()

# ── SECTION 6 — Downstream DE consequences ───────────────────────────────────
st.subheader("📉 Section 6 — How batch effects distort differential expression")

st.markdown("""
Batch effects don't just change how samples look in PCA — they can **inflate or
hide true DE genes** when you compare groups.
""")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown("""
**Without correction / batch modeling:**
- A gene that is high in Batch 1 looks "up in Group A" even if it isn't truly DE
- Batch-driven false positives inflate your DE gene list

**With correction or batch in the model:**
- The batch contribution is subtracted before group comparison
- True biological differences become clearer

**Best practice:**
Include batch as a covariate in your DE model (e.g. in DESeq2 design formula:
`~ batch + groupA`) even if you have already corrected the data for visualisation.
    """)

with col_r:
    # Show a simple 4-sample gene example
    gene_ex = pd.DataFrame({
        "Sample": ["A_batch1", "A_batch2", "B_batch1", "B_batch2"],
        "Group":  ["A", "A", "B", "B"],
        "Batch":  ["Batch1", "Batch2", "Batch1", "Batch2"],
        "Uncorrected": [8.5, 5.5, 6.5, 3.5],
        "Corrected":   [7.0, 7.0, 5.0, 5.0],
    })
    fig_ex = go.Figure()
    for grp, color in [("A", "#3b82f6"), ("B", "#ef4444")]:
        sub = gene_ex[gene_ex["Group"] == grp]
        fig_ex.add_trace(go.Bar(
            name=f"Group {grp} (uncorrected)",
            x=sub["Sample"], y=sub["Uncorrected"],
            marker_color=color, opacity=0.5,
        ))
        fig_ex.add_trace(go.Bar(
            name=f"Group {grp} (corrected)",
            x=sub["Sample"], y=sub["Corrected"],
            marker_color=color, opacity=1.0,
        ))
    fig_ex.update_layout(
        barmode="group",
        title="Example gene: apparent vs true difference",
        yaxis_title="log-CPM",
        height=300,
        margin=dict(t=40),
        legend=dict(font=dict(size=10)),
    )
    st.plotly_chart(fig_ex, use_container_width=True)
    st.caption(
        "The uncorrected values (pale) suggest a large Group A vs B difference. "
        "After removing the batch shift, the true difference is smaller."
    )

st.divider()

# ── SECTION 7 — Cautions ──────────────────────────────────────────────────────
st.subheader("⚠️ Section 7 — Important cautions")

st.error("""
**Never apply batch correction blindly.**

If your batch is confounded with biology (e.g. all controls in Batch 1,
all treated samples in Batch 2), correction will remove real biological signal.
There is no computational fix for a badly designed experiment.
""")

st.warning("""
**Common mistakes to avoid:**

- Correcting when batch and biology are fully confounded
- Applying correction and then forgetting to include batch in the DE model
- Using corrected counts as input to DESeq2 (use raw counts with batch in the design instead)
- Overcorrecting: if batch is weakly present, correction may add noise rather than remove it
""")

st.info("""
**Good practice:**

- Always include batch as a covariate in your DE model design formula
- Use corrected data only for **visualisation** (PCA, heatmaps), not as raw input to DE tools
- Check PCA both before and after correction to confirm biological signal is preserved
- Design studies with balanced batches from the start — correction is a last resort, not a plan
""")

st.divider()

# ── SECTION 8 — Key takeaways ─────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

t1, t2, t3, t4, t5 = st.columns(5)
for col, icon, text in [
    (t1, "🧬", "Batch effects are **technical**, not biological"),
    (t2, "🔭", "PCA reveals batch structure in your data"),
    (t3, "⚖️", "Correction can **improve visualisation** but requires care"),
    (t4, "🎯", "Study design matters more than any correction method"),
    (t5, "📊", "Always model batch in your **DE design formula**"),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:0.9rem;text-align:center;height:120px;display:flex;
align-items:center;justify-content:center;flex-direction:column;">
<div style="font-size:1.6rem">{icon}</div>
<div style="font-size:0.82rem;color:#334155;margin-top:0.3rem">{text}</div>
</div>
""", unsafe_allow_html=True)

st.divider()
col_nav1, col_nav2 = st.columns(2)
with col_nav1:
    st.page_link("pages/04_FDR.py", label="← Lesson 4: FDR", icon="📐")
with col_nav2:
    st.page_link("Home.py", label="Back to course home →", icon="🏠")

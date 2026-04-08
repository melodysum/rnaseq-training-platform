"""
pages/04_FDR.py — Lesson 2: Multiple Testing & FDR
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data, get_paired_columns
from utils.filtering import filter_low_expression
from utils.fdr_demo import (
    simulate_pvalues,
    apply_bh,
    bh_step_table,
    run_paired_de,
)

st.set_page_config(
    page_title="Lesson 4 — FDR",
    page_icon="📐",
    layout="wide",
)

init_session_data()

st.title("📐 Lesson 4 — Multiple Testing & FDR")
st.markdown("""
> **Learning goal:** Understand how Benjamini-Hochberg correction works,
> and — crucially — how your *filtering decision* from Lesson 1 changes the FDR outcome.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔧 Parameters")

    st.subheader("Filtering (linked to Lesson 1)")
    prev = st.session_state.get("filter_results", {})
    min_count = st.slider(
        "Min count per sample",
        0, 100, prev.get("min_count", 10), 1,
    )
    min_samples = st.slider(
        "Min samples passing",
        1, len(counts_raw.columns),
        prev.get("min_samples", max(1, len(counts_raw.columns) // 4)), 1,
    )

    st.subheader("DE significance")
    fdr_cutoff = st.slider("FDR cutoff", 0.01, 0.20, 0.05, 0.01)
    lfc_cutoff = st.slider("Min |log₂FC|", 0.0, 3.0, 1.0, 0.1)

    st.divider()
    st.subheader("BH simulation panel")
    sim_n_genes  = st.slider("Total genes simulated", 500, 20000, 5000, 500)
    sim_n_de     = st.slider("True DE genes in simulation", 0, 500, 100, 10)
    sim_fdr      = st.slider("FDR threshold (simulation)", 0.01, 0.20, 0.05, 0.01)

# ── Run filtering + DE ────────────────────────────────────────────────────────
counts_filtered = filter_low_expression(counts_raw, min_count, min_samples)
n_before = len(counts_raw)
n_after  = len(counts_filtered)
n_removed_filter = n_before - n_after

donors, ctrl_cols, treat_cols = get_paired_columns(metadata)

de_results = None
if ctrl_cols and treat_cols and n_after > 0:
    de_results = run_paired_de(
        counts_filtered, ctrl_cols, treat_cols, fdr_cutoff, lfc_cutoff
    )

# ── Section 1: The multiple testing problem ───────────────────────────────────
st.subheader("🎯 The multiple testing problem")

col_a, col_b = st.columns([3, 2])
with col_a:
    st.markdown(f"""
When you test **{n_after:,} genes** for differential expression, you run
**{n_after:,} statistical tests** — one per gene.

Even if none of them are truly different, a p-value threshold of 0.05 would
*by chance* give you approximately:

> **{int(n_after * 0.05):,} false positives** (5% × {n_after:,} tests)

This is the multiple testing problem. The more tests you run, the more
false positives you expect by chance alone.

**Reducing the number of tested genes (via filtering) directly reduces this problem.**
    """)

with col_b:
    # Show how n_tested changes with filtering
    thresholds = list(range(0, 51, 5))
    n_tested   = [
        len(filter_low_expression(counts_raw, t, min_samples))
        for t in thresholds
    ]
    fig_burden = px.line(
        x=thresholds, y=n_tested,
        labels={"x": "Min count threshold", "y": "Genes tested"},
        title="Testing burden vs filter threshold",
        markers=True,
    )
    fig_burden.add_vline(
        x=min_count, line_dash="dot", line_color="#ef4444",
        annotation_text=f"Current: {min_count}", annotation_position="top right",
    )
    fig_burden.update_layout(height=300, margin=dict(t=40, b=10))
    st.plotly_chart(fig_burden, use_container_width=True)

st.divider()

# ── Section 2: How BH works ───────────────────────────────────────────────────
st.subheader("📊 How Benjamini-Hochberg correction works")

with st.expander("Step-by-step explanation of BH", expanded=True):
    col_exp, col_table = st.columns([2, 3])
    with col_exp:
        st.markdown("""
The **Benjamini-Hochberg (BH)** procedure controls the *False Discovery Rate* —
the expected proportion of significant results that are false.

**The algorithm:**
1. Rank all p-values from smallest (rank 1) to largest (rank m).
2. For each rank k, compute the BH threshold: **(k / m) × α**
3. Find the **largest rank** where p-value ≤ BH threshold.
4. Reject all hypotheses up to that rank.

**Key insight:** The BH threshold gets *stricter* for individual tests as m grows.
With 10,000 genes, a p-value of 0.001 at rank 50 is compared to
(50/10000) × 0.05 = 0.00025 — much stricter than if only 1,000 genes were tested.

This is why **filtering matters**: fewer genes = less strict correction = more power.
        """)
    with col_table:
        # Show BH table using simulated data
        sim_df = simulate_pvalues(sim_n_genes, sim_n_de)
        step_df, critical_rank = bh_step_table(sim_df["pvalue"].values, sim_fdr, show_n=15)
        st.markdown(f"**BH procedure (first 15 rows, m = {sim_n_genes:,}, α = {sim_fdr})**")
        st.dataframe(
            step_df.style
            .format({"p-value": "{:.4f}", "BH threshold (k/m × α)": "{:.4f}"})
            .apply(lambda col: [
                "background-color: #dcfce7" if v else "background-color: #fee2e2"
                for v in col
            ], subset=["Rejected?"])
            .set_properties(**{"font-size": "0.82rem"}),
            use_container_width=True,
            hide_index=True,
        )
        st.caption(f"Critical rank: **{critical_rank}** — all genes up to this rank are rejected.")

st.divider()

# ── Section 3: Interactive BH simulation ─────────────────────────────────────
st.subheader("🧪 Interactive FDR simulation")

sim_df         = simulate_pvalues(sim_n_genes, sim_n_de)
reject, padj   = apply_bh(sim_df["pvalue"].values, sim_fdr)

sim_df["padj"]     = padj
sim_df["rejected"] = reject
sim_df["neg_log10_p"]    = -np.log10(sim_df["pvalue"].clip(lower=1e-300))
sim_df["neg_log10_padj"] = -np.log10(sim_df["padj"].clip(lower=1e-300))

n_sig_sim    = reject.sum()
n_tp_sim     = (reject & sim_df["true_de"]).sum()
n_fp_sim     = (reject & ~sim_df["true_de"]).sum()
fdr_actual   = n_fp_sim / max(n_sig_sim, 1)

m1, m2, m3, m4 = st.columns(4)
m1.metric("Genes tested",    f"{sim_n_genes:,}")
m2.metric("True DE genes",   f"{sim_n_de:,}")
m3.metric("Significant (BH)", f"{n_sig_sim:,}")
m4.metric("Actual FDR",      f"{fdr_actual:.1%}")

col_p, col_padj = st.columns(2)

with col_p:
    fig_p = px.histogram(
        sim_df, x="pvalue", nbins=50,
        color="true_de",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        labels={"pvalue": "Raw p-value", "true_de": "True DE gene"},
        title="Raw p-value distribution",
    )
    fig_p.add_vline(x=sim_fdr, line_dash="dot", line_color="black",
                    annotation_text=f"α = {sim_fdr}")
    fig_p.update_layout(height=320, margin=dict(t=40))
    st.plotly_chart(fig_p, use_container_width=True)

with col_padj:
    fig_padj = px.histogram(
        sim_df, x="padj", nbins=50,
        color="true_de",
        color_discrete_map={True: "#ef4444", False: "#94a3b8"},
        labels={"padj": "BH-adjusted p-value", "true_de": "True DE gene"},
        title="Adjusted p-value distribution",
    )
    fig_padj.add_vline(x=sim_fdr, line_dash="dot", line_color="black",
                       annotation_text=f"FDR = {sim_fdr}")
    fig_padj.update_layout(height=320, margin=dict(t=40))
    st.plotly_chart(fig_padj, use_container_width=True)

st.caption(
    "Red = truly DE genes (ground truth), grey = null genes. "
    "Notice how correction shifts many small raw p-values above the threshold."
)

st.divider()

# ── Section 4: Gene fate — two ways to lose a gene ───────────────────────────
st.subheader("🗺️ Two ways a gene can disappear")

fate_col, fate_viz = st.columns([2, 3])

with fate_col:
    st.markdown(f"""
After all analysis steps, a gene ends up in one of three places:

| Fate | Reason | Stage |
|---|---|---|
| ✅ **Significant DE** | Passed filter + FDR | Both steps |
| 🔴 **Removed by filter** | Low expression | Lesson 1 |
| 🟡 **Tested, not significant** | Failed FDR | Lesson 2 |

With your current settings:
- **{n_removed_filter:,}** genes removed by low-expression filter
- **{n_after:,}** genes entered statistical testing
""")
    if de_results is not None:
        n_sig  = de_results["significant"].sum()
        n_fail = n_after - n_sig
        st.markdown(f"""
- **{n_sig:,}** significant DE genes
- **{n_fail:,}** tested but not significant
        """)

with fate_viz:
    if de_results is not None:
        n_sig  = int(de_results["significant"].sum())
        n_fail = n_after - n_sig

        fig_fate = px.funnel(
            x=[n_before, n_after, n_sig],
            y=["All genes in dataset",
               "After low-expression filter",
               "Significant DE genes"],
            color_discrete_sequence=["#64748b", "#3b82f6", "#22c55e"],
        )
        fig_fate.update_layout(height=320, margin=dict(t=10))
        st.plotly_chart(fig_fate, use_container_width=True)

        # Sankey of fate
        fig_sankey = go.Figure(go.Sankey(
            node=dict(
                label=["All genes", "Removed (filter)",
                       "Tested", "Significant", "Not significant"],
                color=["#94a3b8", "#f87171", "#60a5fa", "#4ade80", "#fbbf24"],
            ),
            link=dict(
                source=[0, 0, 2, 2],
                target=[1, 2, 3, 4],
                value=[n_removed_filter, n_after, n_sig, n_fail],
                color=["#fecaca", "#bfdbfe", "#bbf7d0", "#fde68a"],
            ),
        ))
        fig_sankey.update_layout(height=280, margin=dict(t=10, b=10))
        st.plotly_chart(fig_sankey, use_container_width=True)
        st.caption("Follow genes from the full dataset → filter → statistical testing → significance.")

st.divider()

# ── Section 5: Filtering → FDR link ──────────────────────────────────────────
st.subheader("🔗 How filtering changes your FDR results")

st.markdown("""
The interactive chart below shows how changing the filtering threshold
affects the number of tested genes and (consequently) the number of
significant results.
""")

thresholds_scan = list(range(0, 51, 5))
scan_rows = []
for t in thresholds_scan:
    cf = filter_low_expression(counts_raw, t, min_samples)
    n_t = len(cf)
    n_sig_t = 0
    if ctrl_cols and treat_cols and n_t > 0:
        try:
            de_t = run_paired_de(cf, ctrl_cols, treat_cols, fdr_cutoff, lfc_cutoff)
            n_sig_t = int(de_t["significant"].sum())
        except Exception:
            pass
    scan_rows.append({"Min count": t, "Genes tested": n_t, "Significant DE": n_sig_t})

scan_df = pd.DataFrame(scan_rows)

fig_scan = go.Figure()
fig_scan.add_trace(go.Scatter(
    x=scan_df["Min count"], y=scan_df["Genes tested"],
    name="Genes tested", line=dict(color="#3b82f6"), mode="lines+markers",
))
fig_scan.add_trace(go.Scatter(
    x=scan_df["Min count"], y=scan_df["Significant DE"],
    name="Significant DE", line=dict(color="#22c55e"), mode="lines+markers",
    yaxis="y2",
))
fig_scan.add_vline(x=min_count, line_dash="dot", line_color="#ef4444",
                   annotation_text=f"Current: {min_count}")
fig_scan.update_layout(
    xaxis_title="Min count threshold",
    yaxis=dict(title="Genes tested", color="#3b82f6"),
    yaxis2=dict(title="Significant DE genes", color="#22c55e",
                overlaying="y", side="right"),
    height=380,
    margin=dict(t=20),
    legend=dict(x=0.6, y=0.95),
)
st.plotly_chart(fig_scan, use_container_width=True)
st.caption(
    "As the filter threshold rises, fewer genes are tested (blue), "
    "but the BH correction becomes less severe — often resulting in *more* significant genes (green) up to a point."
)

st.divider()

# ── Volcano of real data ──────────────────────────────────────────────────────
if de_results is not None:
    st.subheader("🌋 Volcano plot — real data DE results")

    de_plot = de_results.reset_index()
    color_map = {
        "Not significant": "#94a3b8",
        "Up":              "#ef4444",
        "Down":            "#3b82f6",
    }
    fig_vol = px.scatter(
        de_plot, x="log2FC", y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={"log2FC": ":.3f", "padj": ":.2e", "direction": False,
                    "neg_log10_padj": False},
        labels={"log2FC": "log₂FC (Treatment / Control)",
                "neg_log10_padj": "−log₁₀(adj. p-value)"},
        category_orders={"direction": ["Up", "Down", "Not significant"]},
        height=480, opacity=0.65,
    )
    fig_vol.update_traces(marker=dict(size=4))
    fig_vol.add_vline(x= lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_vline(x=-lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot",
                      line_color="#64748b", opacity=0.7)
    fig_vol.update_layout(margin=dict(t=10))
    st.plotly_chart(fig_vol, use_container_width=True)

    # Top DE genes table
    st.subheader("📋 Top significant DE genes")
    sig = de_results[de_results["significant"]].sort_values("padj")
    if len(sig) > 0:
        st.dataframe(
            sig.head(30)[["log2FC", "mean_expr", "pvalue", "padj", "direction"]]
            .rename(columns={"log2FC": "log₂FC", "mean_expr": "Mean log-CPM",
                             "pvalue": "p-value", "padj": "adj. p-value",
                             "direction": "Direction"})
            .style.format({"log₂FC": "{:.3f}", "Mean log-CPM": "{:.2f}",
                            "p-value": "{:.2e}", "adj. p-value": "{:.2e}"}),
            use_container_width=True, height=380,
        )
        st.download_button(
            "⬇️ Download full DE results (CSV)",
            data=de_results.reset_index().to_csv(index=False),
            file_name="de_results.csv", mime="text/csv",
        )
    else:
        st.warning("No significant DE genes with current thresholds. Try adjusting the sidebar.")

st.divider()
st.info("""
**Key takeaway from this lesson:**

Filtering and FDR correction are not independent steps — they are directly connected.
Fewer genes tested → less severe BH correction → more power to detect true differences.
The two types of gene loss (filtered out vs FDR failed) are biologically and
statistically different, and you should always report which is which.
""")
st.page_link("Home.py", label="← Back to course home", icon="🏠")

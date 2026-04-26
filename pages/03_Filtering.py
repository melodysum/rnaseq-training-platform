"""
pages/03_Filtering.py — Lesson 1: Low-Expression Filtering
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import (
    filter_low_expression,
    expression_summary,
    check_gene_fate,
    log_cpm,
)

st.set_page_config(
    page_title="Lesson 3 — Filtering",
    page_icon="🔍",
    layout="wide",
)

init_session_data()

# ── Page header ───────────────────────────────────────────────────────────────
st.title("🔍 Lesson 3 — Low-Expression Filtering")
st.markdown("""
> **Learning goal:** Understand *why* we remove low-count genes before analysis,
> and how the threshold choice affects what is kept — and what is thrown away.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
N_GENES    = len(counts_raw)
N_SAMPLES  = len(counts_raw.columns)

# ── Sidebar controls ──────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔧 Filtering Parameters")
    min_count = st.slider(
        "Minimum count per sample",
        0, 100, 10, 1,
        help="A gene needs at least this many reads in enough samples to be kept.",
    )
    min_samples = st.slider(
        "Minimum samples passing threshold",
        1, N_SAMPLES, max(1, N_SAMPLES // 4), 1,
        help="How many samples must exceed the count threshold.",
    )
    st.divider()
    st.caption(
        "💡 Try setting minimum count to 0 (keep everything) "
        "and then raise it — watch how many genes disappear."
    )

# ── Run filtering ─────────────────────────────────────────────────────────────
counts_filtered = filter_low_expression(counts_raw, min_count, min_samples)
n_kept    = len(counts_filtered)
n_removed = N_GENES - n_kept
pct_kept  = n_kept / N_GENES * 100

# Store results in session for FDR page
st.session_state["filter_results"] = {
    "counts_filtered": counts_filtered,
    "min_count":       min_count,
    "min_samples":     min_samples,
    "n_kept":          n_kept,
    "n_removed":       n_removed,
}

# ── Summary metrics ───────────────────────────────────────────────────────────
c1, c2, c3, c4 = st.columns(4)
c1.metric("Total genes",      f"{N_GENES:,}")
c2.metric("Genes retained",   f"{n_kept:,}",    delta=f"{pct_kept:.1f}% kept")
c3.metric("Genes removed",    f"{n_removed:,}", delta=f"-{100-pct_kept:.1f}%", delta_color="inverse")
c4.metric("Samples",          f"{N_SAMPLES}")

st.divider()

# ── Teaching explanation ──────────────────────────────────────────────────────
with st.expander("📖 Why do we filter low-expression genes?", expanded=True):
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("""
**The problem with very low counts**

A gene with counts like 0, 0, 1, 0, 2 across samples isn't giving you reliable
biological signal — it's mostly noise from the sequencing process.

These genes cause two problems:
1. **Inflate variance** — fold changes become huge and meaningless (e.g. 0 → 1 is technically infinite).
2. **Increase testing burden** — every gene tested is another chance for a false positive.
        """)
    with col_b:
        st.markdown("""
**What filtering does**

Filtering removes genes that are too lowly expressed to be reliably measured.

- It is **not arbitrary deletion** — it is a principled preprocessing step.
- The threshold is yours to set: too strict removes real biology,
  too lenient keeps noise.
- A common rule: keep genes with **≥ 10 counts** in at least **as many samples
  as the smallest group size**.

After filtering, every remaining gene has enough data to be tested fairly.
        """)

st.divider()

# ── Distribution plots ────────────────────────────────────────────────────────
st.subheader("📊 Expression distributions: before vs after filtering")

tab1, tab2 = st.tabs(["Histogram", "Per-sample boxplot"])

with tab1:
    rng = np.random.default_rng(42)
    # Sample for speed
    n_plot = min(3000, N_GENES)
    idx_b  = rng.choice(N_GENES,  size=n_plot, replace=False)
    idx_a  = rng.choice(n_kept,   size=min(n_plot, n_kept), replace=False)

    vals_before = np.log2(counts_raw.iloc[idx_b].values.flatten() + 1)
    vals_after  = np.log2(counts_filtered.iloc[idx_a].values.flatten() + 1)

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=vals_before, name="Before filtering",
        opacity=0.55, nbinsx=80, marker_color="#64748b",
    ))
    fig.add_trace(go.Histogram(
        x=vals_after, name="After filtering",
        opacity=0.65, nbinsx=80, marker_color="#3b82f6",
    ))
    fig.update_layout(
        barmode="overlay",
        xaxis_title="log₂(count + 1)",
        yaxis_title="Number of observations",
        height=360,
        margin=dict(t=10),
        legend=dict(x=0.7, y=0.95),
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "The spike of near-zero values (left) represents low/zero-count genes. "
        "Filtering removes this spike, leaving a cleaner distribution."
    )

with tab2:
    # Show first 10 samples for readability
    sample_subset = list(counts_raw.columns[:10])
    lc_before = log_cpm(counts_raw)[sample_subset]
    lc_after  = log_cpm(counts_filtered)[sample_subset] if n_kept > 0 else lc_before

    fig2 = go.Figure()
    for col in sample_subset:
        grp = metadata.loc[col, "groupA"] if col in metadata.index else "unknown"
        color = "#3b82f6" if grp == "control" else "#ef4444"
        fig2.add_trace(go.Box(
            y=lc_before[col], name=col, marker_color=color,
            opacity=0.4, showlegend=False,
        ))
        fig2.add_trace(go.Box(
            y=lc_after[col], name=col + " (filtered)",
            marker_color=color, opacity=0.9, showlegend=False,
        ))
    fig2.update_layout(
        yaxis_title="log₂(CPM + 1)",
        height=380, margin=dict(t=10),
    )
    st.plotly_chart(fig2, use_container_width=True)
    st.caption("Darker = after filtering. Distributions become tighter when noisy low-count genes are removed.")

st.divider()

# ── Gene fate: where did removed genes go? ───────────────────────────────────
st.subheader("🗺️ What happened to removed genes?")

col_l, col_r = st.columns(2)

with col_l:
    st.markdown("""
Genes can leave your analysis at **two different points**:

| Removal stage | Reason | When |
|---|---|---|
| 🔴 **Low-expression filter** | Count too low across samples | Right now, in this step |
| 🟡 **FDR not significant** | Tested but didn't pass correction | Later, in the DE step |

Understanding this distinction matters because:
- A gene lost to filtering was **never tested** statistically.
- A gene lost to FDR **was tested** but didn't show a significant difference.
- These are biologically and statistically very different outcomes.
    """)

with col_r:
    # Pie chart of gene fate
    if n_kept > 0:
        fig_pie = px.pie(
            values=[n_kept, n_removed],
            names=["Retained for testing", "Removed (low expression)"],
            color_discrete_sequence=["#3b82f6", "#f87171"],
            hole=0.45,
        )
        fig_pie.update_layout(height=280, margin=dict(t=10, b=10))
        st.plotly_chart(fig_pie, use_container_width=True)

st.divider()

# ── Gene tracking ─────────────────────────────────────────────────────────────
st.subheader("🔬 Track specific genes")

col_track, col_gt = st.columns(2)

with col_track:
    st.markdown("**Search for any gene**")
    gene_input = st.text_area(
        "Enter gene symbols (one per line or comma-separated)",
        placeholder="e.g.\nHLCS\nNKX2-1\nSNU13",
        height=120,
        key="gene_search",
    )
    if gene_input.strip():
        query_genes = [
            g.strip()
            for g in gene_input.replace(",", "\n").split("\n")
            if g.strip()
        ]
        fate_df = check_gene_fate(counts_raw, counts_filtered, query_genes)
        st.dataframe(
            fate_df.style.apply(
                lambda col: [
                    "background-color: #dcfce7" if "Retained" in v
                    else "background-color: #fee2e2" if "Removed" in v
                    else "background-color: #fef9c3"
                    for v in col
                ],
                subset=["Status"],
            ),
            use_container_width=True,
            hide_index=True,
        )

with col_gt:
    st.markdown("**Ground truth genes** *(known biology)*")
    st.caption(
        "Paste a list of genes you expect to be biologically relevant. "
        "This lets you check whether your filtering threshold is accidentally "
        "removing genes of known interest."
    )
    gt_input = st.text_area(
        "Paste ground truth gene list",
        placeholder="e.g.\nHLCS\nTBC1D13\nH2AZ2\nCLNS1A",
        height=120,
        key="gt_genes",
    )
    if gt_input.strip():
        gt_genes = [
            g.strip()
            for g in gt_input.replace(",", "\n").split("\n")
            if g.strip()
        ]
        gt_fate = check_gene_fate(counts_raw, counts_filtered, gt_genes)

        retained = (gt_fate["Status"].str.contains("Retained")).sum()
        removed  = (gt_fate["Status"].str.contains("Removed")).sum()
        notfound = (gt_fate["Status"].str.contains("Not found")).sum()

        if removed > 0:
            st.warning(
                f"⚠️ {removed} of your ground truth genes are being **removed** "
                f"by the current threshold. Consider lowering the filter."
            )
        if retained > 0:
            st.success(f"✅ {retained} ground truth gene(s) retained.")

        st.dataframe(
            gt_fate.style.apply(
                lambda col: [
                    "background-color: #dcfce7" if "Retained" in v
                    else "background-color: #fee2e2" if "Removed" in v
                    else "background-color: #fef9c3"
                    for v in col
                ],
                subset=["Status"],
            ),
            use_container_width=True,
            hide_index=True,
        )

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — DATA-DRIVEN THRESHOLD (filterByExpr method)
# ══════════════════════════════════════════════════════════════════════════════

from utils.filtering import filter_by_expr, threshold_sweep_retained

st.subheader("📐 Section 5 — Data-Driven Threshold (filterByExpr Method)")
st.markdown("""
The manual threshold above is useful for learning, but in real analysis the
threshold should be derived from the data — not chosen arbitrarily.

**Why arbitrary thresholds are problematic:**  
Setting `min_count = 10` regardless of sequencing depth is inconsistent.  
A gene with 10 reads in a 1M-read library is well-detected;  
the same 10 reads in a 10M-read library represents very low expression.

**filterByExpr logic (edgeR, Chen et al. 2016):**  
The threshold scales with library size:

```
CPM cutoff = min_count / (median_library_size / 1,000,000)
```

A gene is kept if its CPM exceeds this cutoff in at least  
`min_group_size` samples (the smallest experimental group).
""")

with st.expander("⚠️ Why NOT optimise threshold using DE signal?", expanded=False):
    st.markdown("""
A common mistake is to try multiple thresholds, run DE at each, and pick the  
threshold that gives the most significant genes.

**This is circular reasoning:**  
You are using the answer (DE results) to decide the input (which genes to test).  
This inflates false positives because you are effectively tuning the analysis  
to the noise structure of your specific dataset.

**Correct approach:** Choose the threshold using only library sizes and  
group sizes — completely independent of DE results.
""")

col_fe1, col_fe2 = st.columns(2)
with col_fe1:
    fe_min_count = st.slider(
        "filterByExpr min_count",
        1, 50, 10, 1,
        help="Minimum counts in median-sized library. Default 10 is the edgeR recommendation.",
        key="fe_min_count",
    )
    fe_min_total = st.slider(
        "filterByExpr min_total_count",
        5, 50, 15, 1,
        help="Minimum total counts summed across all samples.",
        key="fe_min_total",
    )

counts_fe = filter_by_expr(
    counts_raw, metadata,
    group_col="groupA",
    min_count=fe_min_count,
    min_total_count=fe_min_total,
)

lib_sizes_fe = counts_raw.sum(axis=0)
median_lib_fe = lib_sizes_fe.median()
cpm_cut_fe = fe_min_count / (median_lib_fe / 1e6)

with col_fe2:
    st.metric("Genes retained (filterByExpr)", f"{len(counts_fe):,}")
    st.metric("Genes removed", f"{N_GENES - len(counts_fe):,}")
    st.caption(
        f"Derived CPM cutoff: **{cpm_cut_fe:.3f}** "
        f"(from median library size {median_lib_fe:,.0f})"
    )

st.markdown("#### Threshold sensitivity: genes retained vs threshold")
st.caption(
    "This plot shows how many genes survive each threshold — "
    "**without running DE**. Use it to judge how sensitive the gene count is "
    "to threshold choice, not to maximise a DE signal."
)

sweep_df = threshold_sweep_retained(counts_raw, metadata, group_col="groupA")

fig_sweep = go.Figure()
fig_sweep.add_trace(go.Scatter(
    x=sweep_df["threshold"],
    y=sweep_df["n_retained"],
    mode="lines+markers",
    marker=dict(size=8, color="#2563eb"),
    line=dict(color="#2563eb", width=2),
    text=[f"CPM cut: {r:.3f}<br>Retained: {n:,} ({p}%)"
          for r, n, p in zip(sweep_df["cpm_cutoff"],
                             sweep_df["n_retained"],
                             sweep_df["pct_retained"])],
    hovertemplate="%{text}<extra></extra>",
))
# Mark the current filterByExpr recommendation
rec_row = sweep_df.iloc[(sweep_df["threshold"] - fe_min_count).abs().argsort()[:1]]
fig_sweep.add_trace(go.Scatter(
    x=rec_row["threshold"],
    y=rec_row["n_retained"],
    mode="markers",
    marker=dict(size=14, color="red", symbol="diamond"),
    name=f"Selected threshold ({fe_min_count})",
))
fig_sweep.update_layout(
    xaxis_title="min_count threshold",
    yaxis_title="Genes retained",
    height=350,
    legend=dict(orientation="h", y=-0.25),
)
st.plotly_chart(fig_sweep, use_container_width=True)

st.dataframe(sweep_df, use_container_width=True, hide_index=True)

st.divider()

# ── Key takeaway ──────────────────────────────────────────────────────────────
st.info("""
**Key takeaway from this lesson:**

Filtering is a preprocessing decision, not an analysis step. It determines which genes
will *enter* your statistical test. The threshold you choose here will directly affect
how many genes are tested — and that changes the multiple testing burden, which in turn
changes your FDR results.

👉 **Go to Lesson 4 — FDR** to see exactly how this plays out.
""")

st.page_link("pages/04_FDR.py", label="Continue to Lesson 4 → FDR & Multiple Testing", icon="📐")

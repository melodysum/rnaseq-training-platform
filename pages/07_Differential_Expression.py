"""
pages/07_Differential_Expression.py — Lesson 5: Differential Expression
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data
from utils.filtering import filter_low_expression
from utils.de_analysis import run_de

st.set_page_config(
    page_title="Lesson 7 — Differential Expression",
    page_icon="🧪",
    layout="wide",
)

init_session_data()

st.title("🧪 Lesson 7 — Differential Expression")
st.markdown("""
> **Learning goal:** Understand how gene-level differences between groups are estimated,
> tested, and interpreted in RNA-seq data.
""")

counts_raw = st.session_state["counts"]
metadata   = st.session_state["metadata"]
filter_res = st.session_state.get("filter_results", {})

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔧 DE Parameters")

    # Gene set
    if filter_res:
        counts_use = filter_res["counts_filtered"]
        st.caption(f"Using {len(counts_use):,} filtered genes (Lesson 3).")
    else:
        counts_use = filter_low_expression(counts_raw, 10, 5)
        st.caption(f"Using {len(counts_use):,} genes (default filter).")

    st.divider()

    # Comparison
    groups = sorted(metadata["groupA"].unique())
    ref_group    = st.selectbox("Reference group", groups, index=0)
    target_group = st.selectbox("Target group",    groups,
                                 index=min(1, len(groups)-1))

    st.divider()

    # Optional covariates
    has_batch = "batch" in metadata.columns
    has_donor = "donor" in metadata.columns
    use_batch = st.checkbox("Adjust for batch", value=False,
                             disabled=not has_batch,
                             help="Requires a 'batch' column in metadata.")

    st.divider()

    # Thresholds
    fdr_cutoff = st.slider("FDR cutoff",      0.01, 0.20, 0.05, 0.01)
    lfc_cutoff = st.slider("Min |log₂FC|",    0.0,  3.0,  1.0,  0.1)

    st.caption(
        "⚠️ Educational demo. Uses log-CPM + t-test + BH correction. "
        "Not a substitute for DESeq2 / edgeR / limma-voom."
    )

# ── SECTION 1 — What is DE? ───────────────────────────────────────────────────
st.subheader("🧬 Section 1 — What is differential expression?")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**Differential expression (DE) analysis** asks:

> *"Do genes show systematically different expression levels between two (or more) biological groups?"*

Key points:
- We compare **group means**, not individual samples
- Observed differences might just be **noise** — statistical testing is needed
- Multiple genes are tested simultaneously → we need **FDR correction**
- The result depends on all upstream choices: filtering, normalisation, batch handling
    """)
with col_b:
    st.info("""
**The DE workflow in this lesson:**

1. Filter low-expression genes (Lesson 3)
2. Normalise to log-CPM
3. Optional: adjust for batch effects
4. Gene-wise t-test between groups
5. BH/FDR multiple testing correction
6. Filter by |log₂FC| and FDR threshold
    """)

st.divider()

# ── SECTION 2 — Dataset overview ─────────────────────────────────────────────
st.subheader("📋 Section 2 — Dataset and comparison setup")

c1, c2, c3, c4 = st.columns(4)
c1.metric("Genes tested",  f"{len(counts_use):,}")
c2.metric("Samples",       f"{len(counts_use.columns)}")
c3.metric("Reference",     ref_group)
c4.metric("Target",        target_group)

if ref_group == target_group:
    st.error("Reference and target groups must be different.")
    st.stop()

if has_donor:
    st.info(
        "**Donor column detected.** The test will use a **paired t-test** "
        "(treatment − control per donor), which removes donor-to-donor noise."
    )
elif has_batch and use_batch:
    st.info("**Batch adjustment enabled.** Batch effects will be residualised before testing.")

st.divider()

# ── Run DE ────────────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Running DE analysis…")
def cached_de(counts_key, ref, target, batch, fdr, lfc):
    return run_de(
        counts_use,
        metadata,
        group_col="groupA",
        ref_group=ref,
        target_group=target,
        batch_col="batch" if batch else None,
        donor_col="donor" if has_donor else None,
        fdr_cutoff=fdr,
        lfc_cutoff=lfc,
    )

de_results = cached_de(
    str(counts_use.shape), ref_group, target_group,
    use_batch, fdr_cutoff, lfc_cutoff,
)

# ── SECTION 3 — Results summary ───────────────────────────────────────────────
st.subheader("📊 Section 3 — DE results overview")

n_sig  = de_results["significant"].sum()
n_up   = (de_results["direction"] == "Up").sum()
n_down = (de_results["direction"] == "Down").sum()
test_type = de_results["test_type"].iloc[0]

m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Genes tested",            f"{len(de_results):,}")
m2.metric("Significant DE",          f"{n_sig:,}")
m3.metric(f"⬆ Up in {target_group}", f"{n_up:,}")
m4.metric(f"⬇ Down in {target_group}", f"{n_down:,}")
m5.metric("Test used", test_type)

st.divider()

# ── SECTION 4 — Volcano plot ──────────────────────────────────────────────────
st.subheader("🌋 Section 4 — Volcano plot")

col_vol, col_info = st.columns([3, 1])

with col_vol:
    de_plot = de_results.reset_index()
    color_map = {
        "Not significant": "#94a3b8",
        "Up":   "#ef4444",
        "Down": "#3b82f6",
    }
    y_cap = min(de_plot["neg_log10_padj"].quantile(0.999) * 1.1, 50)

    fig_vol = px.scatter(
        de_plot, x="log2FC", y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={
            "log2FC":          ":.3f",
            "pvalue":          ":.2e",
            "padj":            ":.2e",
            "mean_ref":        ":.2f",
            "mean_target":     ":.2f",
            "direction":       False,
            "neg_log10_padj":  False,
        },
        labels={
            "log2FC":         f"log₂FC ({target_group} / {ref_group})",
            "neg_log10_padj": "−log₁₀(adj. p-value)",
        },
        category_orders={"direction": ["Up", "Down", "Not significant"]},
        height=480,
        opacity=0.65,
    )
    fig_vol.update_traces(marker=dict(size=4))
    fig_vol.add_vline(x= lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_vline(x=-lfc_cutoff, line_dash="dot", line_color="#64748b", opacity=0.7)
    fig_vol.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot",
                      line_color="#64748b", opacity=0.7)
    fig_vol.update_layout(
        yaxis_range=[0, y_cap],
        legend_title="Direction",
        margin=dict(t=10),
    )
    st.plotly_chart(fig_vol, use_container_width=True)

with col_info:
    st.markdown("""
**How to read a volcano plot:**

- **X axis:** log₂ fold change
  - Right = higher in target
  - Left = higher in reference

- **Y axis:** −log₁₀(adj. p-value)
  - Higher = more significant

- **Dashed lines:** your thresholds

- **Top-right / top-left** = biologically and statistically significant

- **Bottom** = not significant regardless of fold change
    """)

st.divider()

# ── SECTION 5 — MA plot ───────────────────────────────────────────────────────
st.subheader("📈 Section 5 — MA plot")

st.markdown("""
The **MA plot** shows fold change vs mean expression.
It reveals whether fold changes are biased at low or high expression levels.
Ideally, the cloud of non-significant genes should centre on log₂FC = 0 across all expression levels.
""")

de_plot["mean_expr"] = (de_plot["mean_ref"] + de_plot["mean_target"]) / 2

fig_ma = px.scatter(
    de_plot, x="mean_expr", y="log2FC",
    color="direction",
    color_discrete_map=color_map,
    hover_name="gene",
    hover_data={"log2FC": ":.3f", "padj": ":.2e", "direction": False},
    labels={
        "mean_expr": f"Mean log-CPM",
        "log2FC":    f"log₂FC ({target_group} / {ref_group})",
    },
    category_orders={"direction": ["Up", "Down", "Not significant"]},
    height=380, opacity=0.55,
)
fig_ma.update_traces(marker=dict(size=3))
fig_ma.add_hline(y=0,           line_dash="solid", line_color="black", opacity=0.3)
fig_ma.add_hline(y= lfc_cutoff, line_dash="dot",   line_color="#64748b", opacity=0.6)
fig_ma.add_hline(y=-lfc_cutoff, line_dash="dot",   line_color="#64748b", opacity=0.6)
fig_ma.update_layout(margin=dict(t=10))
st.plotly_chart(fig_ma, use_container_width=True)
st.caption(
    "Genes at the far left (very low expression) often have unstable fold changes. "
    "This is one reason why low-expression filtering matters."
)

st.divider()

# ── SECTION 6 — Results table ─────────────────────────────────────────────────
st.subheader("📋 Section 6 — Results table")

col_tab, col_search = st.columns([3, 1])

with col_search:
    st.markdown("**Search for a specific gene**")
    gene_query = st.text_input("Gene symbol", placeholder="e.g. HLCS")
    top_n = st.slider("Show top N significant genes", 10, 200, 50, 10)

with col_tab:
    if gene_query.strip():
        query = gene_query.strip()
        if query in de_results.index:
            row = de_results.loc[[query]]
            st.success(f"Found: **{query}**")
            st.dataframe(
                row[["mean_ref", "mean_target", "log2FC", "pvalue", "padj",
                     "significant", "direction"]]
                .rename(columns={
                    "mean_ref":    f"Mean ({ref_group})",
                    "mean_target": f"Mean ({target_group})",
                    "log2FC":      "log₂FC",
                    "pvalue":      "p-value",
                    "padj":        "adj. p-value",
                })
                .style.format({
                    f"Mean ({ref_group})":    "{:.2f}",
                    f"Mean ({target_group})": "{:.2f}",
                    "log₂FC":   "{:.3f}",
                    "p-value":  "{:.2e}",
                    "adj. p-value": "{:.2e}",
                }),
                use_container_width=True,
            )
        else:
            st.warning(f"Gene '{query}' not found in results.")
    else:
        sig = de_results[de_results["significant"]].sort_values("padj")
        if len(sig) == 0:
            st.warning("No significant DE genes. Try adjusting thresholds.")
        else:
            display = sig.head(top_n)[
                ["mean_ref", "mean_target", "log2FC", "pvalue", "padj", "direction"]
            ].rename(columns={
                "mean_ref":    f"Mean ({ref_group})",
                "mean_target": f"Mean ({target_group})",
                "log2FC":      "log₂FC",
                "pvalue":      "p-value",
                "padj":        "adj. p-value",
                "direction":   "Direction",
            })
            st.dataframe(
                display.style.format({
                    f"Mean ({ref_group})":    "{:.2f}",
                    f"Mean ({target_group})": "{:.2f}",
                    "log₂FC":       "{:.3f}",
                    "p-value":      "{:.2e}",
                    "adj. p-value": "{:.2e}",
                }),
                use_container_width=True,
                height=380,
            )

st.download_button(
    "⬇️ Download full DE results (CSV)",
    data=de_results.reset_index().to_csv(index=False),
    file_name="de_results_lesson5.csv",
    mime="text/csv",
)

st.divider()

# ── SECTION 7 — Interpreting DE correctly ─────────────────────────────────────
st.subheader("🎓 Section 7 — Interpreting DE correctly")

col_l, col_r = st.columns(2)
with col_l:
    st.markdown("""
**Statistical significance ≠ biological importance**

| Scenario | Interpretation |
|---|---|
| Small FDR + large log₂FC | Strong candidate |
| Small FDR + tiny log₂FC | Significant but may be negligible biologically |
| Large log₂FC + large FDR | Interesting but unreliable — often low expression |
| Both large | Most compelling |

**log₂FC and FDR answer different questions:**
- **FDR** → "How confident are we this difference is real?"
- **log₂FC** → "How large is the difference?"

You need both to make a good biological interpretation.
    """)

with col_r:
    st.markdown("""
**Upstream decisions affect your DE results**

| Decision | Effect on DE |
|---|---|
| Filtering threshold | Changes the number of tested genes → changes FDR correction |
| Normalisation method | Changes the input to the test |
| Batch correction | Can reduce or inflate apparent differences |
| Paired vs unpaired test | Paired reduces donor noise → more power |

This is why each lesson in this course matters — they all feed into the quality of your final DE results.
    """)

st.divider()

# ── SECTION 8 — Common mistakes ───────────────────────────────────────────────
st.subheader("⚠️ Section 8 — Common mistakes")

st.error("**Using raw p-values without correction.** With thousands of genes tested, raw p < 0.05 gives hundreds of false positives.")
st.warning("**Ignoring low expression.** Fold changes at very low counts are unstable. Filter first.")
st.warning("**Forgetting batch effects.** Unmodelled batch can inflate DE gene lists dramatically.")
st.info("**Using corrected counts as input to DESeq2.** Batch-corrected values are for visualisation — use raw counts with batch in the design formula instead.")
st.info("**Treating the volcano plot as the final answer.** It is a summary, not a conclusion. Follow up top hits with functional context.")

st.divider()

# ── SECTION 9 — Takeaways ─────────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

cols = st.columns(5)
msgs = [
    ("🧬", "Groups not samples", "DE compares group means, not individual measurements."),
    ("🔢", "Both metrics matter", "FDR and log₂FC together make a result interpretable."),
    ("📊", "Volcano = summary", "Plots help exploration — they are not statistical proof."),
    ("🔗", "Upstream choices matter", "Filtering, normalisation, and batch all shape your DE results."),
    ("🎯", "Design is everything", "Paired designs, good metadata, and balanced batches improve power."),
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
st.page_link("pages/05_Exploratory_Analysis_PCA.py",
             label="← Lesson 5: Exploratory Analysis", icon="📊")
st.page_link("Home.py", label="Back to course home →", icon="🏠")

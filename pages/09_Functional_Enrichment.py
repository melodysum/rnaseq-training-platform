"""
pages/09_Functional_Enrichment.py
Lesson 9 — Functional Enrichment (GO / GSEA)

Educational enrichment analysis using built-in toy pathways.
Not a substitute for clusterProfiler, fgsea, or MSigDB-based pipelines.
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

from utils.data_loader import init_session_data, load_demo_data
from utils.enrichment_utils import (
    TOY_GENE_SETS, run_ora, run_gsea_like, get_running_sum, demo_ranked_genes,
)

st.set_page_config(
    page_title="Lesson 9 — Functional Enrichment",
    page_icon="🧩",
    layout="wide",
)
init_session_data()

st.title("🧩 Lesson 9 — Functional Enrichment (GO / GSEA)")
st.caption(
    "Educational implementation using built-in toy pathways. "
    "Not a substitute for clusterProfiler, fgsea, or MSigDB-based analysis."
)

# ════════════════════════════════════════════════════════════════════════════════
# 1. CONCEPT EXPLANATION
# ════════════════════════════════════════════════════════════════════════════════

with st.expander("📖 What is functional enrichment analysis?", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
### Over-Representation Analysis (ORA)
ORA asks: *are genes from a biological pathway appearing in my DE gene list
more often than expected by chance?*

**How it works:**
1. Take your list of significant DE genes.
2. For each pathway in a database (e.g. GO, KEGG), count how many pathway genes
   appear in your list.
3. Use Fisher's exact test to ask whether that overlap is larger than expected
   given the background gene universe.
4. Correct for multiple testing (Benjamini-Hochberg FDR).

**Strength:** Simple, fast, interpretable.  
**Limitation:** Treats all significant genes equally — ignores fold-change magnitude.
The choice of background universe significantly affects results.
        """)
    with col2:
        st.markdown("""
### Gene Set Enrichment Analysis (GSEA)
GSEA asks: *does a pathway tend to appear at the top (or bottom) of a ranked
gene list, even if individual genes don't pass a significance cutoff?*

**How it works:**
1. Rank **all** genes by a statistic (e.g. log2FC, signed p-value).
2. Walk down the ranked list; when a gene belongs to the pathway, step up;
   otherwise step down.
3. The enrichment score (ES) is the maximum deviation of this running sum.
4. Pathways with genes concentrated at the top are positively enriched;
   at the bottom, negatively enriched.

**Strength:** Uses the full gene list — more sensitive than ORA.  
**Limitation:** More complex to interpret; permutation-based significance is
computationally expensive (educational version uses simpler scoring).
        """)
    st.info("""
**Why does the gene universe matter?**  
In ORA, the background (universe) defines what "expected by chance" means.
Using only the genes you measured (not all 20,000 human genes) gives more
realistic enrichment statistics. A narrow universe inflates, while an
overly broad one deflates enrichment p-values.
    """)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 2. DATA SOURCE
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("📂 Data Source")

counts   = st.session_state["counts"]
metadata = st.session_state["metadata"]
source   = st.session_state.get("data_source", "demo")
de_res   = st.session_state.get("de_results", None)

# Check whether usable DE results exist in session
have_de = (
    de_res is not None
    and isinstance(de_res, pd.DataFrame)
    and "log2FC" in de_res.columns
    and len(de_res) > 0
)

if have_de:
    st.success(
        f"✅ Found DE results in session ({len(de_res):,} genes). "
        "You can use them as the gene list source below.",
        icon="🧪",
    )
else:
    st.info(
        "No DE results found in session state. "
        "A ranked gene list will be generated from the "
        f"{'built-in demo' if source == 'demo' else 'uploaded'} data "
        "using log2FC between the two groups.",
        icon="📊",
    )
    # Generate demo ranking from current data
    try:
        de_res = demo_ranked_genes(counts, metadata)
        have_de = True
        st.caption(
            f"Demo ranking: {de_res['group2'].iloc[0]} vs {de_res['group1'].iloc[0]} "
            f"({len(de_res):,} genes ranked by log2FC)."
        )
    except Exception as e:
        st.error(f"Could not generate demo ranking: {e}")
        have_de = False

st.caption(
    f"**Pathway library:** {len(TOY_GENE_SETS)} built-in educational toy pathways "
    f"({sum(len(v) for v in TOY_GENE_SETS.values())} total genes across all sets). "
    "These are for teaching only — not a real GO/MSigDB database."
)

st.divider()

if not have_de:
    st.error("Cannot run enrichment analysis without a gene list. "
             "Please run Lesson 7 — Differential Expression first.")
    st.stop()

# ════════════════════════════════════════════════════════════════════════════════
# 3. GENE LIST SOURCE SELECTION
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🎯 Choose Gene List Source")

gene_list_options = []
if have_de and "log2FC" in de_res.columns:
    gene_list_options += [
        "Significant DE genes (padj < 0.05, |log2FC| ≥ 1)",
        "Top 100 upregulated genes",
        "Top 100 downregulated genes",
    ]
gene_list_options.append("Manual gene list (enter below)")

source_choice = st.radio("Gene list source:", gene_list_options, horizontal=True)

if source_choice == "Significant DE genes (padj < 0.05, |log2FC| ≥ 1)":
    if "padj" in de_res.columns and "significant" in de_res.columns:
        sig = de_res[de_res["significant"]].copy()
    elif "padj" in de_res.columns:
        sig = de_res[(de_res["padj"] < 0.05) & (de_res["log2FC"].abs() >= 1)].copy()
    else:
        sig = de_res[de_res["log2FC"].abs() >= 1].head(200).copy()
    gene_list_for_ora = sig.index.tolist() if sig.index.name == "gene" else sig["gene"].tolist() if "gene" in sig.columns else sig.index.tolist()
    st.caption(f"{len(gene_list_for_ora)} significant DE genes selected.")

elif source_choice == "Top 100 upregulated genes":
    top = de_res.sort_values("log2FC", ascending=False).head(100)
    gene_list_for_ora = top.index.tolist() if top.index.name == "gene" else top["gene"].tolist() if "gene" in top.columns else top.index.tolist()
    st.caption(f"Top {len(gene_list_for_ora)} upregulated genes by log2FC.")

elif source_choice == "Top 100 downregulated genes":
    top = de_res.sort_values("log2FC", ascending=True).head(100)
    gene_list_for_ora = top.index.tolist() if top.index.name == "gene" else top["gene"].tolist() if "gene" in top.columns else top.index.tolist()
    st.caption(f"Top {len(gene_list_for_ora)} downregulated genes by log2FC.")

else:
    manual_input = st.text_area(
        "Enter gene symbols (one per line or comma-separated):",
        placeholder="STAT1\nMX1\nISG15\nIFIT1",
        height=120,
    )
    gene_list_for_ora = [
        g.strip().upper()
        for g in manual_input.replace(",", "\n").splitlines()
        if g.strip()
    ]
    st.caption(f"{len(gene_list_for_ora)} genes entered.")

if len(gene_list_for_ora) < 3:
    st.warning("Please select or enter at least 3 genes to run enrichment analysis.")
    st.stop()

# Build ranked list for GSEA-like (always from full DE results)
if "gene" in de_res.columns:
    ranked_series = de_res.set_index("gene")["log2FC"].sort_values(ascending=False)
elif de_res.index.name == "gene" or de_res.index.dtype == object:
    ranked_series = de_res["log2FC"].sort_values(ascending=False)
else:
    ranked_series = de_res["log2FC"].sort_values(ascending=False)

lfc_weights = ranked_series.abs()

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 4. ORA
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("📊 Over-Representation Analysis (ORA)")
st.caption(
    "Fisher's exact test + Benjamini-Hochberg FDR correction. "
    "Background universe = all genes measured in the current dataset."
)

universe = counts.index.tolist()

with st.spinner("Running ORA…"):
    ora_res = run_ora(gene_list_for_ora, TOY_GENE_SETS, universe=universe)

if ora_res.empty:
    st.warning("No pathway overlaps found. Try a larger gene list or different source.")
else:
    # Display table
    display_cols = ["pathway", "overlap", "pathway_size", "query_size",
                    "gene_ratio", "pvalue", "padj", "overlap_genes"]
    display_cols = [c for c in display_cols if c in ora_res.columns]
    st.dataframe(
        ora_res[display_cols].style.format({
            "gene_ratio": "{:.3f}",
            "pvalue":     "{:.4f}",
            "padj":       "{:.4f}",
        }),
        use_container_width=True,
        height=280,
    )

    # Bar plot — top pathways by -log10(padj)
    plot_df = ora_res.copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-10))
    plot_df["label"] = plot_df["pathway"].str.replace("_", " ").str.title()
    plot_df = plot_df.sort_values("neg_log10_padj", ascending=True).tail(10)

    fig_ora = px.bar(
        plot_df,
        x="neg_log10_padj",
        y="label",
        orientation="h",
        color="gene_ratio",
        color_continuous_scale="Blues",
        labels={
            "neg_log10_padj": "-log₁₀(adjusted p-value)",
            "label":          "Pathway",
            "gene_ratio":     "Gene Ratio",
        },
        title="Top Enriched Pathways (ORA) — Educational Toy Pathways",
    )
    fig_ora.add_vline(x=-np.log10(0.05), line_dash="dash", line_color="red",
                      annotation_text="padj = 0.05")
    fig_ora.update_layout(height=400)
    st.plotly_chart(fig_ora, use_container_width=True)
    st.caption(
        "⚠️ These are educational toy pathways, not a real GO or KEGG database. "
        "Results are illustrative only."
    )

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 5. GSEA-LIKE RUNNING SUM ENRICHMENT
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("📈 Ranked Enrichment Analysis (GSEA-like)")
st.caption(
    "Genes ranked by log2FC. Running-sum enrichment score walks down the ranked list. "
    "Labeled as **educational GSEA-like implementation** — not GSEA permutation testing."
)

with st.spinner("Running GSEA-like scoring…"):
    gsea_res = run_gsea_like(ranked_series.index, TOY_GENE_SETS, lfc_weights)

if gsea_res.empty:
    st.warning("No enrichment results produced.")
else:
    display_cols_g = ["pathway", "ES", "direction", "n_hits", "pathway_size"]
    display_cols_g = [c for c in display_cols_g if c in gsea_res.columns]

    gsea_display = gsea_res[display_cols_g].copy()
    gsea_display["pathway"] = gsea_display["pathway"].str.replace("_", " ").str.title()

    st.dataframe(
        gsea_display.style.format({"ES": "{:.4f}"}),
        use_container_width=True,
        height=300,
    )

    # Enrichment plot for user-selected pathway
    st.markdown("#### Enrichment Plot — Select a Pathway")
    st.markdown(
        "The enrichment plot shows the running sum as we walk down the ranked gene list. "
        "A peak at the top indicates the pathway genes are concentrated among "
        "the most upregulated genes."
    )

    pathway_options = gsea_res["pathway"].tolist()
    pathway_labels  = [p.replace("_", " ").title() for p in pathway_options]
    selected_label  = st.selectbox(
        "Select pathway to plot:", pathway_labels, index=0
    )
    selected_pathway = pathway_options[pathway_labels.index(selected_label)]
    pathway_genes    = TOY_GENE_SETS[selected_pathway]

    rs = get_running_sum(ranked_series.index, pathway_genes, lfc_weights)
    ranks = list(range(1, len(rs) + 1))

    # Build plotly figure
    fig_es = go.Figure()

    # Running sum line
    fig_es.add_trace(go.Scatter(
        x=ranks, y=rs,
        mode="lines",
        name="Running Enrichment Score",
        line=dict(color="#2563eb", width=2),
    ))

    # Zero line
    fig_es.add_hline(y=0, line_dash="dash", line_color="gray", line_width=1)

    # Tick marks where pathway genes hit
    hit_positions = [i + 1 for i, g in enumerate(ranked_series.index)
                     if g in set(pathway_genes)]
    for pos in hit_positions:
        fig_es.add_vline(x=pos, line_color="rgba(220,38,38,0.25)", line_width=0.8)

    # Mark ES peak
    rs_array = np.array(rs)
    if abs(rs_array.max()) >= abs(rs_array.min()):
        peak_idx = int(np.argmax(rs_array))
    else:
        peak_idx = int(np.argmin(rs_array))
    fig_es.add_trace(go.Scatter(
        x=[ranks[peak_idx]], y=[rs_array[peak_idx]],
        mode="markers",
        name=f"ES = {rs_array[peak_idx]:.3f}",
        marker=dict(color="red", size=10, symbol="diamond"),
    ))

    fig_es.update_layout(
        title=(
            f"Enrichment Plot: {selected_label}<br>"
            "<sup>Educational GSEA-like implementation — Simulated data</sup>"
        ),
        xaxis_title="Gene Rank (log2FC descending)",
        yaxis_title="Running Enrichment Score",
        height=420,
        legend=dict(orientation="h", y=-0.2),
    )
    st.plotly_chart(fig_es, use_container_width=True)

    es_val = gsea_res[gsea_res["pathway"] == selected_pathway]["ES"].values[0]
    direction = gsea_res[gsea_res["pathway"] == selected_pathway]["direction"].values[0]
    n_hits_val = gsea_res[gsea_res["pathway"] == selected_pathway]["n_hits"].values[0]

    st.info(
        f"**{selected_label}** — ES = {es_val:.4f} | Direction: {direction} | "
        f"Pathway genes found in ranked list: {n_hits_val}/{len(pathway_genes)}\n\n"
        "Vertical red lines show where pathway genes land in the ranked list. "
        "A concentrated cluster at the top or bottom produces a clear ES peak."
    )

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 6. INTERPRETATION NOTES
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("💡 Interpreting Enrichment Results")
st.markdown("""
**ORA results:**
- A low adjusted p-value suggests the pathway is over-represented in your gene list.
- Gene Ratio (overlap / pathway size) indicates what fraction of pathway members
  are in your DE gene list — a high ratio with a small overlap can still be
  meaningful for small pathways.
- Always check the overlap genes — sometimes a handful of highly connected hub
  genes drive enrichment across many pathways.

**GSEA-like results:**
- A positive ES means pathway genes cluster near the top of the ranked list
  (associated with upregulation in the target group).
- A negative ES means pathway genes cluster near the bottom (downregulated).
- The enrichment plot shape matters: a sharp peak close to rank 1 is stronger
  evidence than a shallow peak mid-list.

**General cautions:**
- These toy pathways are for education only. Real analysis requires GO, KEGG,
  Reactome, or MSigDB databases.
- Enrichment is highly sensitive to the gene universe and ranking statistic.
- Biological plausibility should guide interpretation, not just p-values.
""")

st.divider()

# ══════════════════════════════════════════════════════════════════════════════
# ADVANCED: GSEA WITH PERMUTATION TESTING
# ══════════════════════════════════════════════════════════════════════════════

from utils.enrichment_utils import run_gsea_permutation, rank_by_statistic
from utils.de_analysis import run_de

st.subheader("🧬 Advanced: GSEA with Permutation Testing & NES")

with st.expander("📖 How this differs from the GSEA-like analysis above", expanded=True):
    st.markdown("""
| | GSEA-like (above) | Permutation GSEA (here) |
|---|---|---|
| **Ranking** | log2FC | Test statistic (t-stat) |
| **Null distribution** | None | Gene-label permutation |
| **Score** | ES only | ES + **NES** (normalised) |
| **p-value** | None | Permutation-based |
| **FDR** | None | Benjamini-Hochberg |
| **Comparable across pathways** | ❌ | ✅ |

**Why use test statistic for ranking instead of log2FC?**  
The t-statistic integrates both effect size and variance.  
A gene with logFC = 3 in 3 samples has a lower t-stat than a gene with  
logFC = 2 in 20 samples — the latter is much more reliable.  
log2FC alone inflates the ranking of noisy, under-replicated genes.

**Why do ORA and GSEA give different results?**  
ORA throws away all sub-threshold genes. If a pathway has 12 genes all  
slightly below FDR 0.05, ORA finds nothing. GSEA sees that all 12 genes  
are near the top of the ranked list and correctly identifies the pathway.
""")

with st.spinner("Running differential expression for GSEA ranking…"):
    de_for_gsea = run_de(counts, metadata, fdr_cutoff=0.05, lfc_cutoff=0.5)

# Add a 'stat' column (t-statistic approximation) for ranking
if "stat" not in de_for_gsea.columns:
    de_for_gsea["stat"] = de_for_gsea["log2FC"] / (
        de_for_gsea["log2FC"].abs().mean() /
        (-np.log10(de_for_gsea["pvalue"].clip(1e-300))).clip(lower=0.01)
    ).clip(lower=0.01)

ranked_for_perm = rank_by_statistic(
    de_for_gsea,
    stat_col="stat",
    lfc_col="log2FC",
    pval_col="pvalue",
    method="signed_log",
)

col_pg1, col_pg2 = st.columns([2, 1])
with col_pg1:
    n_perm = st.select_slider(
        "Number of permutations",
        options=[100, 500, 1000],
        value=500,
        help="More permutations = more stable p-values. 1000 is standard.",
    )
with col_pg2:
    min_set = st.number_input("Min pathway size", 3, 20, 5, key="min_set_perm")

run_perm = st.button("▶ Run Permutation GSEA", type="primary")

if run_perm or "perm_gsea_result" in st.session_state:
    if run_perm:
        with st.spinner(f"Running {n_perm} permutations…"):
            perm_result = run_gsea_permutation(
                ranked_for_perm,
                gene_sets=TOY_GENE_SETS,
                n_permutations=n_perm,
                random_state=42,
                min_set_size=min_set,
            )
            st.session_state["perm_gsea_result"] = perm_result
    else:
        perm_result = st.session_state["perm_gsea_result"]

    if perm_result.empty:
        st.warning("No pathways passed the minimum size filter.")
    else:
        display_perm = perm_result[[
            "pathway", "n_genes_in_list", "ES", "NES",
            "pvalue", "padj", "direction", "leading_edge_genes"
        ]].copy()
        display_perm["pathway"] = display_perm["pathway"].str.replace("_", " ").str.title()
        display_perm["significant"] = display_perm["padj"] < 0.05

        def _style_sig(row):
            if row["padj"] < 0.05:
                return ["background-color: #dcfce7"] * len(row)
            return [""] * len(row)

        st.dataframe(
            display_perm.style
                .apply(_style_sig, axis=1)
                .format({"ES": "{:.4f}", "NES": "{:.4f}",
                         "pvalue": "{:.4f}", "padj": "{:.4f}"}),
            use_container_width=True,
            height=350,
            hide_index=True,
        )

        sig_perm = (perm_result["padj"] < 0.05).sum()
        st.info(
            f"**{sig_perm}** pathway(s) with FDR < 0.05 after permutation testing. "
            "Green rows are significant. NES > 0 = pathway upregulated; "
            "NES < 0 = pathway downregulated."
        )

        st.markdown("#### NES bar chart")
        fig_nes = go.Figure()
        colors = ["#16a34a" if v > 0 else "#dc2626" for v in perm_result["NES"]]
        pway_labels = [p.replace("_", " ").title() for p in perm_result["pathway"]]
        fig_nes.add_trace(go.Bar(
            y=pway_labels,
            x=perm_result["NES"],
            orientation="h",
            marker_color=colors,
            text=[f"padj={p:.3f}" for p in perm_result["padj"]],
            textposition="outside",
        ))
        fig_nes.add_vline(x=0, line_color="black", line_width=1)
        fig_nes.update_layout(
            xaxis_title="NES (Normalized Enrichment Score)",
            yaxis=dict(autorange="reversed"),
            height=max(300, len(perm_result) * 35),
        )
        st.plotly_chart(fig_nes, use_container_width=True)
else:
    st.info("Click **▶ Run Permutation GSEA** to compute NES and permutation p-values.")

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# NAVIGATION
# ════════════════════════════════════════════════════════════════════════════════

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/08_Clustering_Heatmaps.py",
                 label="← Lesson 8: Clustering & Heatmaps", icon="🗺️")
with col_n2:
    st.page_link("pages/10_Public_Data.py",
                 label="Lesson 10: Public Data & Reproducibility →", icon="🌐")

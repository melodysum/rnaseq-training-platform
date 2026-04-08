"""
pages/01_RNAseq_Foundations.py — Lesson 1: RNA-seq Foundations & Experimental Design
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

st.set_page_config(
    page_title="Lesson 1 — RNA-seq Foundations",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 Lesson 1 — RNA-seq Foundations & Experimental Design")
st.markdown("""
> **Learning goal:** Understand where RNA-seq data come from, what core sequencing
> concepts mean, and why good experimental design matters *before* any analysis begins.
""")

# ── SECTION 1 — Why foundations matter ───────────────────────────────────────
st.subheader("🏗️ Section 1 — Why foundations matter")

c1, c2, c3 = st.columns(3)
for col, icon, title, body in [
    (c1, "📐", "Design first",
     "RNA-seq analysis does not begin with statistics — it begins with how samples were designed, collected, and sequenced."),
    (c2, "🔗", "Everything is connected",
     "Misunderstanding FASTQ quality, replicates, or depth leads to confusion later in filtering, PCA, and DE analysis."),
    (c3, "⚠️", "Damage control is limited",
     "Bad experimental design decisions often cannot be fully repaired computationally, no matter how sophisticated the downstream method."),
]:
    col.markdown(f"""
<div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;
padding:1rem;height:140px;">
<div style="font-size:1.4rem">{icon}</div>
<div style="font-weight:600;margin:0.3rem 0">{title}</div>
<div style="font-size:0.85rem;color:#475569">{body}</div>
</div>
""", unsafe_allow_html=True)

st.divider()

# ── SECTION 2 — What is RNA-seq? ─────────────────────────────────────────────
st.subheader("🧬 Section 2 — What is RNA-seq?")

st.markdown("""
**RNA-seq** (RNA sequencing) measures the abundance of RNA transcripts in a biological sample.
It gives a genome-wide snapshot of which genes are active and at what level.

- **Bulk RNA-seq** measures *average* expression across thousands or millions of cells
- Reads are later aligned to a reference and counted to produce an **expression matrix**
- The expression matrix (genes × samples) is the starting point for all downstream analysis
""")

# Workflow diagram
st.markdown("**The RNA-seq workflow:**")
stages = [
    ("🧫", "Biological\nSample", "Tissue, cells, or organoids"),
    ("🔬", "Library\nPreparation", "RNA extraction, fragmentation, cDNA synthesis, adapter ligation"),
    ("💻", "Sequencing", "Short reads generated (FASTQ files)"),
    ("📁", "FASTQ\nFiles", "Raw sequencing data with quality scores"),
    ("🗺️", "Quantification", "Reads aligned and counted per gene/transcript"),
    ("📊", "Count\nMatrix", "Genes × samples table of raw counts"),
    ("🔍", "Downstream\nAnalysis", "Filtering, normalisation, DE, enrichment"),
]

fig_flow = go.Figure()
n = len(stages)
for i, (icon, label, tooltip) in enumerate(stages):
    x = i / (n - 1)
    fig_flow.add_trace(go.Scatter(
        x=[x], y=[0.5],
        mode="markers+text",
        marker=dict(size=40, color="#3b82f6", symbol="circle"),
        text=[icon],
        textposition="middle center",
        textfont=dict(size=18),
        hovertext=tooltip,
        hoverinfo="text",
        showlegend=False,
    ))
    fig_flow.add_annotation(
        x=x, y=0.18, text=label.replace("\n", "<br>"),
        showarrow=False, font=dict(size=10, color="#334155"),
        align="center",
    )
    if i < n - 1:
        fig_flow.add_annotation(
            x=(i + 0.5) / (n - 1), y=0.5,
            text="→", showarrow=False,
            font=dict(size=20, color="#94a3b8"),
        )

fig_flow.update_layout(
    height=180, margin=dict(l=10, r=10, t=10, b=40),
    xaxis=dict(visible=False, range=[-0.05, 1.05]),
    yaxis=dict(visible=False, range=[0, 1]),
    plot_bgcolor="white",
)
st.plotly_chart(fig_flow, use_container_width=True)
st.caption("Hover over each stage for more detail.")

st.divider()

# ── SECTION 3 — FASTQ anatomy ─────────────────────────────────────────────────
st.subheader("📄 Section 3 — FASTQ file anatomy")

col_fastq, col_qual = st.columns([3, 2])

with col_fastq:
    st.markdown("""
Every sequencing read is stored as a **4-line FASTQ record**:
    """)
    st.code("""@SRR123456.1 read_1 length=150          ← Read identifier
ATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATG  ← Nucleotide sequence
+                                                  ← Separator (always +)
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!        ← Quality string (Phred)
""", language="text")

    with st.expander("🔍 What do the quality characters mean?"):
        st.markdown("""
Each character encodes a **Phred quality score** (Q):

| Q score | Error probability | ASCII character (example) |
|---|---|---|
| Q10 | 1 in 10 | + |
| Q20 | 1 in 100 | 5 |
| Q30 | 1 in 1,000 | ? |
| Q40 | 1 in 10,000 | I |

- Quality is encoded as ASCII character: `Q = ASCII_value - 33`
- `!` = Q0 (worst), `I` = Q40 (excellent)
- Quality typically **drops at the end of reads** due to cumulative sequencing chemistry errors
        """)

with col_qual:
    st.markdown("**Quality score explorer**")
    q_score = st.slider("Phred Q score", min_value=0, max_value=40, value=30)
    error_prob = 10 ** (-q_score / 10)
    accuracy   = (1 - error_prob) * 100

    st.metric("Error probability", f"1 in {int(1/error_prob):,}" if error_prob > 0 else "< 1 in 10,000")
    st.metric("Base call accuracy", f"{accuracy:.4f}%")

    q_vals = list(range(0, 41))
    probs  = [10 ** (-q / 10) for q in q_vals]
    fig_q = px.line(
        x=q_vals, y=probs,
        labels={"x": "Q score", "y": "Error probability"},
        title="Q score → error probability",
        log_y=True,
    )
    fig_q.add_vline(x=q_score, line_dash="dot", line_color="#ef4444",
                    annotation_text=f"Q{q_score}", annotation_position="top right")
    fig_q.update_layout(height=260, margin=dict(t=40, b=10))
    st.plotly_chart(fig_q, use_container_width=True)

st.divider()

# ── SECTION 4 — Core sequencing concepts ─────────────────────────────────────
st.subheader("📚 Section 4 — Core sequencing concepts")

concepts = {
    "Read length": {
        "definition": "The number of base pairs sequenced per read (e.g. 75 bp, 150 bp).",
        "why": "Longer reads improve alignment confidence and splice-junction detection.",
        "misconception": "Longer is not always better — 150 bp is often sufficient for gene-level DE analysis.",
    },
    "Sequencing depth": {
        "definition": "Total number of reads generated per sample (e.g. 20M, 50M reads).",
        "why": "Greater depth improves detection of lowly expressed genes and reduces sampling noise.",
        "misconception": "More depth cannot replace biological replication — deeper sequencing of 2 samples still gives you n=2.",
    },
    "Library": {
        "definition": "A collection of DNA fragments prepared from one sample, ready for sequencing.",
        "why": "Library quality (RNA integrity, adapter ligation) directly affects data quality.",
        "misconception": "A technical replicate (two libraries from the same RNA) is not the same as a biological replicate.",
    },
    "Single-end vs paired-end": {
        "definition": "Single-end reads one end of each fragment; paired-end reads both ends.",
        "why": "Paired-end improves alignment and is preferred for transcript isoform analysis.",
        "misconception": "Paired-end does not double biological replication — it just gives two reads per fragment.",
    },
    "Biological replicates": {
        "definition": "Independent biological samples from the same condition (e.g. 4 separate mice from the same group).",
        "why": "Essential for estimating within-group variance and for any statistical test to be valid.",
        "misconception": "Having n=1 or n=2 per group is not sufficient for reliable DE — variability cannot be estimated.",
    },
    "Technical replicates": {
        "definition": "Repeated measurements of the same biological sample (e.g. sequencing the same RNA twice).",
        "why": "Useful for assessing platform reproducibility, but should not substitute biological replicates.",
        "misconception": "Technical replicates cannot tell you about biological variability between individuals.",
    },
}

tabs = st.tabs(list(concepts.keys()))
for tab, (concept, info) in zip(tabs, concepts.items()):
    with tab:
        col_def, col_misc = tab.columns(2)
        with col_def:
            st.markdown(f"**Definition:** {info['definition']}")
            st.markdown(f"**Why it matters:** {info['why']}")
        with col_misc:
            st.warning(f"⚠️ **Common misconception:**\n\n{info['misconception']}")

st.divider()

# ── SECTION 5 — Experimental design simulator ────────────────────────────────
st.subheader("🧪 Section 5 — Experimental design simulator")

st.markdown("Explore how design choices affect your study's power and cost.")

col_ctrl, col_output = st.columns([1, 2])

with col_ctrl:
    n_groups     = st.slider("Number of biological groups",  2, 6, 2)
    n_reps       = st.slider("Biological replicates per group", 2, 10, 4)
    reads_per_m  = st.slider("Sequencing depth per sample (M reads)", 5, 100, 20)
    cost_per_m   = st.number_input("Estimated cost per million reads (£)", 1.0, 20.0, 4.0, 0.5)
    has_batch    = st.checkbox("Add batch structure (2 batches)")

n_samples     = n_groups * n_reps
total_reads_m = n_samples * reads_per_m
total_cost    = total_reads_m * cost_per_m

if has_batch:
    samples_b1 = n_samples // 2
    samples_b2 = n_samples - samples_b1
    batch_balanced = (samples_b1 == samples_b2) and (n_reps % 2 == 0)

with col_output:
    m1, m2, m3 = st.columns(3)
    m1.metric("Total samples",    f"{n_samples}")
    m2.metric("Total reads (M)",  f"{total_reads_m:,}")
    m3.metric("Est. cost",        f"£{total_cost:,.0f}")

    # Design interpretation
    issues = []
    if n_reps < 3:
        issues.append("⚠️ **Too few replicates.** n < 3 per group severely limits power and variance estimation.")
    if reads_per_m < 10:
        issues.append("⚠️ **Low depth.** < 10M reads may miss lowly expressed genes.")
    if reads_per_m > 60 and n_reps < 4:
        issues.append("⚠️ **High depth but low replication.** Consider redistributing sequencing budget to more replicates.")
    if has_batch and not batch_balanced:
        issues.append("⚠️ **Unbalanced batch structure.** Some groups may be confounded with batch.")
    if not issues:
        st.success("✅ **Solid design.** Good replication, reasonable depth, and balanced structure.")
    for issue in issues:
        st.warning(issue)

    # Visual: replicates vs depth trade-off
    rep_range   = list(range(2, 11))
    depth_range = [max(5, int(total_reads_m / (n_groups * r))) for r in rep_range]
    fig_td = go.Figure()
    fig_td.add_trace(go.Scatter(
        x=rep_range, y=depth_range,
        mode="lines+markers", name="Depth per sample",
        line=dict(color="#3b82f6"),
        marker=dict(size=8),
    ))
    fig_td.add_vline(x=n_reps, line_dash="dot", line_color="#ef4444",
                     annotation_text=f"Current: {n_reps} reps", annotation_position="top right")
    fig_td.add_hline(y=20, line_dash="dot", line_color="#22c55e",
                     annotation_text="20M recommended", annotation_position="right")
    fig_td.update_layout(
        xaxis_title="Replicates per group",
        yaxis_title="Reads per sample (M)",
        title="Replicates vs depth trade-off (fixed budget)",
        height=300, margin=dict(t=40),
    )
    st.plotly_chart(fig_td, use_container_width=True)

st.divider()

# ── SECTION 6 — Replication vs depth ─────────────────────────────────────────
st.subheader("⚖️ Section 6 — Replication vs depth")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**Design A — High depth, few replicates**
- 2 replicates per group
- 50M reads per sample
- Total cost: similar to Design B

❌ Problems:
- Cannot reliably estimate within-group variability
- Statistical tests have very low power
- One outlier sample ruins the comparison
- Cannot tell if a difference is real or just noise
    """)
    st.error("**Not recommended for DE analysis**")

with col_b:
    st.markdown("""
**Design B — Balanced replication**
- 4–5 replicates per group
- 15–20M reads per sample
- Same or lower total cost

✅ Benefits:
- Robust variance estimation per gene
- Higher statistical power for DE
- More resistant to sample outliers
- Suitable for most standard DE workflows
    """)
    st.success("**Preferred design for most RNA-seq studies**")

st.info("""
**Key insight:** Adding one more biological replicate usually improves your DE analysis
more than doubling the sequencing depth. Depth matters for detection of rare transcripts,
but replication is essential for statistical inference.
""")

st.divider()

# ── SECTION 7 — Common design failures ───────────────────────────────────────
st.subheader("⚠️ Section 7 — Common design failures")

failures = [
    ("🚫", "No biological replication",
     "A single sample per group makes it statistically impossible to estimate variance or run a valid test.",
     "error"),
    ("🚫", "No proper control group",
     "Without an appropriate untreated or baseline control, fold changes have no reference point.",
     "error"),
    ("⚠️", "Batch confounded with condition",
     "All treated samples in Batch 1, all controls in Batch 2 — computational correction cannot untangle these.",
     "warning"),
    ("⚠️", "Poor or missing metadata",
     "If you cannot identify which sample is which after sequencing, the data becomes uninterpretable.",
     "warning"),
    ("ℹ️", "Relying only on RNA-seq to validate RNA-seq",
     "Significant DE results should ideally be validated by orthogonal methods (qPCR, protein).",
     "info"),
    ("ℹ️", "Treating technical replicates as biological",
     "Inflates the apparent n and gives false confidence. These are not independent observations.",
     "info"),
]

for icon, title, body, level in failures:
    getattr(st, level)(f"**{icon} {title}**\n\n{body}")

st.divider()

# ── SECTION 8 — Key takeaways ─────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

cols = st.columns(5)
msgs = [
    ("🧬", "RNA-seq measures transcripts", "It gives a genome-wide snapshot of gene activity."),
    ("📄", "FASTQ = sequence + quality", "Quality scores guide filtering and trimming decisions."),
    ("👥", "Replicates are essential", "Without biological replication, statistical inference is not valid."),
    ("⚖️", "Replication > depth", "More replicates at moderate depth beats fewer replicates at high depth."),
    ("🏗️", "Design cannot be fixed later", "Good experimental design is more powerful than any computational correction."),
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
col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("Home.py", label="← Back to course home", icon="🏠")
with col_n2:
    st.page_link("pages/02_Quantification_Import_Annotation.py",
                 label="Lesson 2: Quantification & Annotation →", icon="📦")

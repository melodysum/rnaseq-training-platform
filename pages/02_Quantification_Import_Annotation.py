"""
pages/02_Quantification_Import_Annotation.py — Lesson 2: Quantification, Import & Annotation
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from utils.data_loader import init_session_data

st.set_page_config(
    page_title="Lesson 2 — Quantification & Annotation",
    page_icon="📦",
    layout="wide",
)

init_session_data()

st.title("📦 Lesson 2 — Quantification, Import & Annotation")
st.markdown("""
> **Learning goal:** Understand how sequencing reads become expression matrices,
> and why metadata and annotation are essential before any downstream analysis.
""")

# ── SECTION 1 — From FASTQ to count matrix ────────────────────────────────────
st.subheader("🗺️ Section 1 — From FASTQ to expression matrix")

stages = [
    ("📁", "FASTQ", "Raw reads + quality scores from the sequencer"),
    ("📚", "Reference\ntranscriptome", "Known transcript sequences for the organism"),
    ("🔢", "Quantification", "Map reads to transcripts; estimate abundance"),
    ("🧮", "Transcript\nabundances", "Per-transcript counts, TPM, effective length"),
    ("🔗", "tx2gene\nsummarisation", "Aggregate transcripts → gene-level values"),
    ("📊", "Gene-level\nmatrix", "Genes × samples count matrix — ready for analysis"),
]

fig = go.Figure()
n = len(stages)
for i, (icon, label, tip) in enumerate(stages):
    x = i / (n - 1)
    fig.add_trace(go.Scatter(
        x=[x], y=[0.55], mode="markers+text",
        marker=dict(size=44, color="#6366f1", symbol="circle"),
        text=[icon], textposition="middle center",
        textfont=dict(size=18),
        hovertext=tip, hoverinfo="text", showlegend=False,
    ))
    fig.add_annotation(
        x=x, y=0.15, text=label.replace("\n", "<br>"),
        showarrow=False, font=dict(size=10, color="#334155"), align="center",
    )
    if i < n - 1:
        fig.add_annotation(
            x=(i + 0.5) / (n - 1), y=0.55, text="→",
            showarrow=False, font=dict(size=20, color="#94a3b8"),
        )
fig.update_layout(
    height=170, margin=dict(l=10, r=10, t=10, b=50),
    xaxis=dict(visible=False, range=[-0.05, 1.05]),
    yaxis=dict(visible=False, range=[0, 1]),
    plot_bgcolor="white",
)
st.plotly_chart(fig, use_container_width=True)
st.caption("Hover over each stage for details.")

st.divider()

# ── SECTION 2 — Alignment vs pseudoalignment ──────────────────────────────────
st.subheader("⚡ Section 2 — Alignment vs pseudoalignment")

col_al, col_ps = st.columns(2)

with col_al:
    st.markdown("""
**Traditional alignment** *(e.g. STAR, HISAT2)*

Maps every read base-by-base to a reference genome.

✅ Precise positional information  
✅ Detects novel splice junctions  
✅ Genome-level output (useful for variant calling)  
❌ Slow and memory-intensive  
❌ Requires a full genome reference  
    """)

with col_ps:
    st.markdown("""
**Pseudoalignment** *(e.g. kallisto, salmon)*

Identifies which transcripts are *compatible* with each read — no base-level mapping needed.

✅ Very fast (10–100× faster)  
✅ Low memory requirements  
✅ Transcript-level estimates built in  
❌ No precise genomic position  
❌ Requires a transcriptome reference  
    """)

comparison_df = pd.DataFrame({
    "Feature":        ["Speed", "Memory", "Output level", "Novel splicing", "Best for"],
    "STAR/HISAT2":    ["Slow", "High (30+ GB)", "Genome", "Yes", "Full characterisation"],
    "kallisto/salmon":["Very fast", "Low (< 4 GB)", "Transcript", "No", "Quantification-focused DE"],
})
st.dataframe(comparison_df, use_container_width=True, hide_index=True)

st.divider()

# ── SECTION 3 — Counts, TPM, effective length ────────────────────────────────
st.subheader("🔢 Section 3 — Counts, TPM, and effective length")

col_def, col_interactive = st.columns([2, 3])

with col_def:
    st.markdown("""
**Counts** — number of reads assigned to a gene/transcript.
- Integer values, reflect raw read support
- Used directly by DESeq2/edgeR (which model count distributions)
- *Do not normalise for gene length or library size*

**TPM** (Transcripts Per Million)
- Normalises for both gene length AND library size
- Useful for comparing expression levels within or across samples
- *Not recommended as input to count-based DE tools*

**Effective length**
- Adjusts for the fragment size distribution
- A shorter effective length means fewer positions where a fragment could originate
- Affects abundance estimation: `estimated counts = TPM × effective_length / 1e6 × lib_size`
    """)

with col_interactive:
    st.markdown("**Interactive: how length and depth affect TPM**")
    n_transcripts = 3
    demo_data = pd.DataFrame({
        "Transcript": ["TxA", "TxB", "TxC"],
        "Raw counts": [500, 200, 50],
        "Length (bp)": [1000, 500, 250],
    })

    col_s1, col_s2 = st.columns(2)
    with col_s1:
        tx_len_A = st.slider("Length of TxA (bp)", 200, 3000, 1000, 100)
    with col_s2:
        tx_len_B = st.slider("Length of TxB (bp)", 200, 3000, 500, 100)

    frag_len = st.slider("Mean fragment length (bp)", 50, 400, 180, 10)

    lengths = [max(1, tx_len_A - frag_len + 1),
               max(1, tx_len_B - frag_len + 1),
               max(1, 250 - frag_len + 1)]
    counts  = [500, 200, 50]
    rpk     = [c / (l / 1000) for c, l in zip(counts, lengths)]
    tpm     = [r / sum(rpk) * 1e6 for r in rpk]

    result_df = pd.DataFrame({
        "Transcript":       ["TxA", "TxB", "TxC"],
        "Raw counts":       counts,
        "Eff. length (bp)": [int(l) for l in lengths],
        "TPM":              [f"{t:.1f}" for t in tpm],
    })
    st.dataframe(result_df, use_container_width=True, hide_index=True)
    st.caption(
        "Adjust transcript lengths and fragment length to see how TPM changes "
        "even when raw counts stay the same."
    )

st.divider()

# ── SECTION 4 — Import checklist ─────────────────────────────────────────────
st.subheader("✅ Section 4 — Importing quantified data")

col_check, col_example = st.columns([2, 3])

with col_check:
    st.markdown("""
Before analysis, your data import should satisfy:
    """)

    checklist = [
        ("✅", "Quantification files present for all samples"),
        ("✅", "Sample names in files match metadata"),
        ("✅", "Metadata file includes groupA and batch columns"),
        ("✅", "Reference transcriptome/annotation version recorded"),
        ("✅", "File paths are structured and reproducible"),
        ("⚠️", "If sample names differ, analysis will fail silently or give wrong results"),
    ]
    for icon, item in checklist:
        st.markdown(f"{icon} {item}")

with col_example:
    st.markdown("**Example: name mismatch breaks analysis**")

    ok_df = pd.DataFrame({
        "counts.csv column": ["D01_control", "D01_treatment", "D02_control"],
        "metadata sample_name": ["D01_control", "D01_treatment", "D02_control"],
        "Match?": ["✅", "✅", "✅"],
    })
    bad_df = pd.DataFrame({
        "counts.csv column": ["sample1", "sample2", "sample3"],
        "metadata sample_name": ["D01_control", "D01_treatment", "D02_control"],
        "Match?": ["❌", "❌", "❌"],
    })

    st.markdown("**✅ Good — names match:**")
    st.dataframe(ok_df, use_container_width=True, hide_index=True)
    st.markdown("**❌ Bad — names mismatch:**")
    st.dataframe(bad_df, use_container_width=True, hide_index=True)
    st.error("Name mismatches cause metadata to detach from samples — your analysis loses biological meaning.")

st.divider()

# ── SECTION 5 — Transcript-to-gene summarisation ─────────────────────────────
st.subheader("🔗 Section 5 — Transcript-to-gene summarisation")

st.markdown("""
One gene often has **multiple transcripts** (splice variants).
For standard DE analysis, we usually work at **gene level** — we need to aggregate.
""")

col_tx, col_gene = st.columns(2)

with col_tx:
    st.markdown("**Transcript-level (from kallisto/salmon):**")
    tx_df = pd.DataFrame({
        "Transcript ID":  ["ENST001.1", "ENST001.2", "ENST001.3"],
        "Gene ID":        ["ENSG001"] * 3,
        "Gene symbol":    ["BRCA1"] * 3,
        "TPM":            [120.3, 45.2, 8.1],
        "Est. counts":    [540, 210, 38],
    })
    st.dataframe(tx_df, use_container_width=True, hide_index=True)

with col_gene:
    st.markdown("**After tximport / gene-level aggregation:**")
    gene_df = pd.DataFrame({
        "Gene symbol": ["BRCA1"],
        "Gene ID":     ["ENSG001"],
        "Sum counts":  [788],
        "Sum TPM":     [173.6],
    })
    st.dataframe(gene_df, use_container_width=True, hide_index=True)
    st.info(
        "tximport aggregates transcript-level values to gene level, "
        "accounting for effective length during the summarisation."
    )

st.divider()

# ── SECTION 6 — Annotation ────────────────────────────────────────────────────
st.subheader("🏷️ Section 6 — Gene annotation and IDs")

col_ann, col_map = st.columns([2, 3])

with col_ann:
    st.markdown("""
**Layers of gene identity:**

| ID type | Example | Notes |
|---|---|---|
| Transcript ID | ENST00000367770.8 | Specific isoform + version |
| Gene ID | ENSG00000139618.19 | Stable Ensembl gene + version |
| Gene symbol | BRCA2 | Human-readable but not always stable |

**Why annotation matters:**
- Version mismatches cause IDs to not match between files
- Gene symbols can map to multiple IDs, or change over time
- Non-model organisms may have incomplete or no annotation
    """)

with col_map:
    st.markdown("**Annotation mapping example:**")
    map_df = pd.DataFrame({
        "Transcript ID":     ["ENST00000367770.8", "ENST00000544455.6"],
        "→ Gene ID":         ["ENSG00000139618.19", "ENSG00000139618.19"],
        "→ Gene symbol":     ["BRCA2", "BRCA2"],
        "→ Chromosome":      ["13", "13"],
        "Annotation source": ["Ensembl v109", "Ensembl v109"],
    })
    st.dataframe(map_df, use_container_width=True, hide_index=True)

    st.warning("""
**Annotation pitfalls:**
- Gene symbols are NOT always unique — some symbols map to multiple genes
- Always record the Ensembl version you used
- Strip version suffixes (`.8`, `.19`) before matching IDs across databases
    """)

st.divider()

# ── SECTION 7 — Metadata ──────────────────────────────────────────────────────
st.subheader("📋 Section 7 — Metadata as the backbone of analysis")

metadata = st.session_state.get("metadata")

if metadata is not None and len(metadata) > 0:
    st.success("✅ Your uploaded metadata is shown below.")
    st.dataframe(metadata.head(10), use_container_width=True)
else:
    st.info("No data uploaded — showing example metadata structure.")
    example_meta = pd.DataFrame({
        "sample_name":   ["D01_control", "D01_treatment", "D02_control", "D02_treatment"],
        "groupA":        ["control", "treatment", "control", "treatment"],
        "donor":         ["D01", "D01", "D02", "D02"],
        "batch":         ["batch1", "batch1", "batch2", "batch2"],
        "sex":           ["F", "F", "M", "M"],
        "age":           [37, 37, 40, 40],
    })
    st.dataframe(example_meta, use_container_width=True, hide_index=True)

st.markdown("""
| Column | Used for | Lesson |
|---|---|---|
| `groupA` | Group comparison / DE | Lessons 3, 4, 7 |
| `batch` | Batch modelling / correction | Lessons 4, 6 |
| `donor` | Paired design | Lessons 4, 7 |
| `sex`, `age` | Covariates / interpretation | Lessons 5, 7 |
""")

st.divider()

# ── SECTION 8 — Common pitfalls ───────────────────────────────────────────────
st.subheader("⚠️ Section 8 — Common pitfalls")

st.error("**Sample name mismatch.** Counts and metadata silently disconnect if names differ.")
st.error("**Using TPM as DE input.** DESeq2 and edgeR require raw counts — TPM removes the count distribution properties these tools depend on.")
st.warning("**Ignoring transcript aggregation.** Running DE on transcript-level without aggregation mixes isoforms inappropriately.")
st.warning("**Ignoring annotation versions.** IDs from different Ensembl versions may not match.")
st.info("**Assuming gene symbols are unique.** Some symbols map to paralogues or pseudogenes.")

st.divider()

# ── SECTION 9 — Takeaways ─────────────────────────────────────────────────────
st.subheader("📌 Key takeaways")

cols = st.columns(5)
msgs = [
    ("🗺️", "Quantification bridges reads and analysis", "FASTQ → count matrix requires several careful steps."),
    ("🔢", "Counts ≠ TPM", "They answer different questions and have different use cases."),
    ("🔗", "Aggregate transcripts to genes", "Gene-level analysis is standard for DE workflows."),
    ("🏷️", "Annotation version matters", "Always record which reference and version you used."),
    ("📋", "Metadata is everything", "Good metadata enables every downstream interpretation."),
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
    st.page_link("pages/01_RNAseq_Foundations.py",
                 label="← Lesson 1: RNA-seq Foundations", icon="🔬")
with col_n2:
    st.page_link("pages/03_Filtering.py",
                 label="Lesson 3: Filtering →", icon="🔍")

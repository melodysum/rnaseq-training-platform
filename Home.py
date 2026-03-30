"""
Home.py — Landing page for the RNA-seq Interactive Training Platform.
Run with: streamlit run Home.py
"""

import streamlit as st
import pandas as pd
from utils.data_loader import validate_and_parse, init_session_data, load_demo_data

st.set_page_config(
    page_title="RNA-seq Training Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Global CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
    html, body, [class*="css"] { font-family: 'Inter', sans-serif; }

    .hero-title {
        font-size: 2.6rem;
        font-weight: 700;
        color: #0f172a;
        line-height: 1.2;
        margin-bottom: 0.4rem;
    }
    .hero-sub {
        font-size: 1.1rem;
        color: #475569;
        max-width: 680px;
        line-height: 1.6;
    }
    .lesson-card {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 12px;
        padding: 1.2rem 1.4rem;
        margin-bottom: 0.8rem;
        transition: box-shadow 0.2s;
    }
    .lesson-card:hover { box-shadow: 0 4px 16px rgba(0,0,0,0.07); }
    .lesson-title { font-size: 1.05rem; font-weight: 600; color: #1e293b; }
    .lesson-desc  { font-size: 0.88rem; color: #64748b; margin-top: 0.2rem; }
    .badge-available {
        background: #dcfce7; color: #166534;
        border-radius: 999px; padding: 2px 10px;
        font-size: 0.75rem; font-weight: 600;
    }
    .badge-coming {
        background: #f1f5f9; color: #94a3b8;
        border-radius: 999px; padding: 2px 10px;
        font-size: 0.75rem; font-weight: 600;
    }
    .path-step {
        display: flex; align-items: flex-start;
        margin-bottom: 0.7rem;
    }
    .step-num {
        background: #3b82f6; color: white;
        border-radius: 50%; width: 28px; height: 28px;
        display: flex; align-items: center; justify-content: center;
        font-size: 0.8rem; font-weight: 700;
        flex-shrink: 0; margin-right: 0.7rem; margin-top: 2px;
    }
    .step-text { font-size: 0.92rem; color: #334155; }
    .divider { border: none; border-top: 1px solid #e2e8f0; margin: 1.5rem 0; }
    .section-label {
        font-size: 0.75rem; font-weight: 700;
        text-transform: uppercase; letter-spacing: 0.08em;
        color: #94a3b8; margin-bottom: 0.6rem;
    }
</style>
""", unsafe_allow_html=True)

# ── Initialise session data ───────────────────────────────────────────────────
init_session_data()

# ── Hero ──────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="hero-title">🧬 RNA-seq Interactive Training Platform</div>
<div class="hero-sub">
    An interactive, step-by-step learning environment for understanding
    transcriptomic data analysis — from raw counts to biological interpretation.
</div>
""", unsafe_allow_html=True)

st.markdown('<hr class="divider">', unsafe_allow_html=True)

# ── Two-column layout ─────────────────────────────────────────────────────────
col_left, col_right = st.columns([3, 2], gap="large")

with col_left:

    # ── Available lessons ─────────────────────────────────────────────────────
    st.markdown('<div class="section-label">Course Content</div>', unsafe_allow_html=True)

    lessons = [
        {
            "icon": "🔬",
            "title": "Lesson 1 — RNA-seq Foundations & Experimental Design",
            "desc": "Understand how RNA-seq data are generated, what key sequencing concepts mean, and why good experimental design matters before downstream analysis.",
            "page": "pages/01_RNAseq_Foundations.py",
            "available": True,
        },
        {
            "icon": "📦",
            "title": "Lesson 2 — Quantification, Import & Annotation",
            "desc": "Learn how RNA-seq reads are converted into transcript and gene-level expression data, and how metadata and annotation make downstream analysis possible.",
            "page": "pages/02_Quantification_Import_Annotation.py",
            "available": True,
        },
        {
            "icon": "🔍",
            "title": "Lesson 3 — Low-Expression Filtering",
            "desc": "Understand why low-count genes are removed before analysis, and how the filtering threshold affects what genes are retained.",
            "page": "pages/03_Filtering.py",
            "available": True,
        },
        {
            "icon": "📐",
            "title": "Lesson 4 — Multiple Testing & FDR",
            "desc": "Learn how Benjamini-Hochberg correction works, why filtering changes your FDR results, and how to distinguish statistical from biological gene loss.",
            "page": "pages/04_FDR.py",
            "available": True,
        },
        {
            "icon": "📊",
            "title": "Lesson 5 — Exploratory Analysis & PCA",
            "desc": "Visualise sample relationships, identify outliers, and interpret principal component analysis in RNA-seq data.",
            "page": "pages/05_Exploratory_Analysis_PCA.py",
            "available": True,
        },
        {
            "icon": "🔄",
            "title": "Lesson 6 — Batch Correction",
            "desc": "Explore how technical batch effects arise, how to detect them with PCA, and how to correct them before differential expression analysis.",
            "page": "pages/06_Batch_Correction.py",
            "available": True,
        },
        {
            "icon": "🧪",
            "title": "Lesson 7 — Differential Expression",
            "desc": "Run a proper DE analysis, understand the statistical model, and interpret results correctly.",
            "page": "pages/07_Differential_Expression.py",
            "available": True,
        },
        {
            "icon": "🗺️",
            "title": "Lesson 8 — Clustering & Heatmaps",
            "desc": "Learn how to cluster genes and samples, build expression heatmaps, and identify gene modules from DE results.",
            "page": "pages/08_Clustering_Heatmaps.py",
            "available": False,
        },
        {
            "icon": "🧩",
            "title": "Lesson 9 — Functional Enrichment (GO / GSEA)",
            "desc": "Run GO over-representation and GSEA analyses to interpret the biological meaning of DE gene lists.",
            "page": "pages/09_Functional_Enrichment.py",
            "available": False,
        },
        {
            "icon": "🌐",
            "title": "Lesson 10 — Public Data & Reproducible Analysis",
            "desc": "Download data from GEO, build reproducible workflows, and document your analysis for sharing and publication.",
            "page": "pages/10_Public_Data.py",
            "available": False,
        },
    ]

    for lesson in lessons:
        badge = (
            '<span class="badge-available">Available</span>'
            if lesson["available"]
            else '<span class="badge-coming">Coming soon</span>'
        )
        st.markdown(f"""
        <div class="lesson-card">
            <div class="lesson-title">{lesson["icon"]} {lesson["title"]} &nbsp; {badge}</div>
            <div class="lesson-desc">{lesson["desc"]}</div>
        </div>
        """, unsafe_allow_html=True)
        if lesson["page"]:
            label = f"Open {lesson['title']} →" if lesson["available"] else "Preview (coming soon)"
            st.page_link(lesson["page"], label=label)

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    # ── Recommended path ──────────────────────────────────────────────────────
    st.markdown('<div class="section-label">Recommended Learning Path</div>', unsafe_allow_html=True)
    steps = [
        ("Start here", "Read <b>Lesson 1 — Foundations</b> to understand what RNA-seq measures and how experimental design shapes everything downstream."),
        ("How data are made", "Study <b>Lesson 2 — Quantification</b> to understand how FASTQ reads become a count matrix."),
        ("Clean your data", "Open <b>Lesson 3 — Filtering</b> and explore how threshold choices affect retained genes."),
        ("Understand FDR", "Move to <b>Lesson 4 — FDR</b> to see how filtering and multiple testing are connected."),
        ("Explore structure", "Use <b>Lesson 5 — Exploratory Analysis</b> to inspect PCA and detect batch effects or outliers."),
        ("Correct batch", "Open <b>Lesson 6 — Batch Correction</b> if technical variation is present."),
        ("Run DE", "Complete <b>Lesson 7 — Differential Expression</b> to compare groups and interpret results."),
    ]
    for i, (label, text) in enumerate(steps, 1):
        st.markdown(f"""
        <div class="path-step">
            <div class="step-num">{i}</div>
            <div class="step-text"><b>{label}:</b> {text}</div>
        </div>
        """, unsafe_allow_html=True)

with col_right:

    # ── Data upload ───────────────────────────────────────────────────────────
    st.markdown('<div class="section-label">Your Data</div>', unsafe_allow_html=True)

    source = st.session_state.get("data_source", "demo")
    if source == "demo":
        st.info(
            "**Currently using built-in demo data.**\n\n"
            "Upload your own files below to analyse your dataset.",
            icon="📂",
        )
    else:
        st.success("✅ Using your uploaded data.", icon="✅")
        if st.button("Switch back to demo data"):
            counts, metadata = load_demo_data()
            st.session_state["counts"]   = counts
            st.session_state["metadata"] = metadata
            st.session_state["data_source"] = "demo"
            st.session_state.pop("filter_results", None)
            st.rerun()

    with st.expander("📁 Upload your own data", expanded=(source == "demo")):
        st.markdown("""
**counts.csv** — raw count matrix
```
gene_symbol, D01_control, D01_treatment, ...
IFNA6,       25,          40,            ...
CAD,          0,           9,            ...
```

**metadata.csv** — one row per sample
```
sample_name,   groupA,    donor, batch,  sex, age
D01_control,   control,   D01,   batch1, F,   37
D01_treatment, treatment, D01,   batch1, F,   37
```
        """)

        counts_file   = st.file_uploader("Upload counts.csv",   type="csv", key="upload_counts")
        metadata_file = st.file_uploader("Upload metadata.csv", type="csv", key="upload_meta")

        if counts_file and metadata_file:
            counts, metadata, err = validate_and_parse(counts_file, metadata_file)
            if err:
                st.error(f"❌ {err}")
            else:
                st.session_state["counts"]   = counts
                st.session_state["metadata"] = metadata
                st.session_state["data_source"] = "uploaded"
                st.session_state.pop("filter_results", None)
                st.success(
                    f"✅ Loaded **{len(counts):,} genes** × **{len(counts.columns)} samples**."
                )

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    # ── Dataset summary ───────────────────────────────────────────────────────
    st.markdown('<div class="section-label">Current Dataset</div>', unsafe_allow_html=True)
    counts   = st.session_state["counts"]
    metadata = st.session_state["metadata"]

    groups = metadata["groupA"].value_counts()
    st.markdown(f"""
| | |
|---|---|
| **Genes** | {len(counts):,} |
| **Samples** | {len(counts.columns)} |
| **Groups** | {', '.join(f'{g} ({n})' for g, n in groups.items())} |
| **Source** | {'Demo data' if source == 'demo' else 'Uploaded'} |
""")

    if "donor" in metadata.columns:
        st.caption(f"Paired design detected: {metadata['donor'].nunique()} donors.")
    if "batch" in metadata.columns:
        st.caption(f"Batch info available: {metadata['batch'].nunique()} batches.")

    st.markdown('<hr class="divider">', unsafe_allow_html=True)

    # ── Quick navigation ──────────────────────────────────────────────────────
    st.markdown('<div class="section-label">Jump to a lesson</div>', unsafe_allow_html=True)
    st.page_link("pages/01_RNAseq_Foundations.py",           label="🔬 Lesson 1 — Foundations",           icon="▶")
    st.page_link("pages/02_Quantification_Import_Annotation.py", label="📦 Lesson 2 — Quantification",        icon="▶")
    st.page_link("pages/03_Filtering.py",                     label="🔍 Lesson 3 — Filtering",              icon="▶")
    st.page_link("pages/04_FDR.py",                           label="📐 Lesson 4 — FDR",                    icon="▶")
    st.page_link("pages/05_Exploratory_Analysis_PCA.py",      label="📊 Lesson 5 — Exploratory Analysis",  icon="▶")
    st.page_link("pages/06_Batch_Correction.py",              label="🔄 Lesson 6 — Batch Correction",      icon="▶")
    st.page_link("pages/07_Differential_Expression.py",       label="🧪 Lesson 7 — Differential Expression", icon="▶")

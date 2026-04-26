"""
pages/10_Public_Data.py
Lesson 10 — Public Data & Reproducible Analysis
"""

import streamlit as st
import pandas as pd
from utils.data_loader import init_session_data

st.set_page_config(
    page_title="Lesson 10 — Public Data & Reproducibility",
    page_icon="🌐",
    layout="wide",
)
init_session_data()

st.title("🌐 Lesson 10 — Public Data & Reproducible Analysis")
st.caption(
    "How public RNA-seq datasets are structured, how to choose the right files, "
    "and how to make your analysis reproducible."
)

# ════════════════════════════════════════════════════════════════════════════════
# 1. GEO / ARRAYEXPRESS CONCEPTS
# ════════════════════════════════════════════════════════════════════════════════

with st.expander("📖 What are GEO and ArrayExpress?", expanded=True):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
### NCBI GEO (Gene Expression Omnibus)
GEO is the primary public repository for high-throughput gene expression data,
maintained by NCBI. It accepts data from microarrays, RNA-seq, ChIP-seq, and
many other genomics platforms.

Every submission receives a **GSE accession** (e.g. `GSE148036`) for the study,
and **GSM accessions** for individual samples. The **GPL accession** describes
the platform (e.g. Illumina HiSeq 2500).

When a paper says "data are available at GEO under accession GSE…", this is
where you look. You can search at [ncbi.nlm.nih.gov/geo](https://www.ncbi.nlm.nih.gov/geo/).
        """)
    with col2:
        st.markdown("""
### EBI ArrayExpress / BioStudies
ArrayExpress (now integrated into EBI BioStudies) is the European counterpart
to GEO. Many studies deposit data in both repositories. Accessions follow the
format `E-MTAB-XXXX`.

**Sequence Read Archive (SRA):** raw FASTQ reads submitted to GEO are stored
in SRA and can be downloaded with `prefetch` / `fastq-dump` from the SRA Toolkit.
European data is mirrored in EBI's ENA (European Nucleotide Archive).

Both repositories are free, open-access, and searchable by organism, tissue,
disease, or keyword.
        """)

    st.markdown("---")
    st.markdown("""
### What files are typically available?

| File type | What it contains | When to use it |
|-----------|-----------------|----------------|
| **Raw FASTQ** | Sequencing reads; available via SRA/ENA | Re-align from scratch; most rigorous |
| **Count matrix** | Gene × sample read counts; in supplementary files | Re-analysis from counts; faster start |
| **Processed / normalised** | TPM, FPKM, RPKM tables; in supplementary files | Visualisation, exploration (not DE) |
| **Series matrix** | GEO-formatted metadata + expression summary | Microarray; also contains sample metadata |
| **SOFT / MINiML** | GEO metadata in structured format | Automated metadata parsing |

For bulk RNA-seq **re-analysis**, the gene-level count matrix in supplementary files
is usually the best starting point — it avoids re-running the full alignment
pipeline while preserving the integer count values needed for DESeq2/edgeR.
    """)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 2. CASE STUDY EXPLORER
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🔬 Real Public Study Case Studies")
st.caption(
    "Two real bulk RNA-seq GEO studies summarised for teaching purposes. "
    "Details are stored locally — no internet access required."
)

CASE_STUDIES = [
    {
        "accession": "GSE148036",
        "title": "Transcriptomic profiling of human alveolar macrophages during Mycobacterium tuberculosis infection",
        "organism": "Homo sapiens",
        "study_question": (
            "How do human alveolar macrophages respond transcriptionally to "
            "Mycobacterium tuberculosis infection, and what host pathways are "
            "dysregulated during active versus latent TB?"
        ),
        "experimental_design": (
            "Bulk RNA-seq of primary human alveolar macrophages isolated from "
            "healthy donors (n = 10) and patients with active pulmonary TB (n = 10). "
            "Two conditions per donor: uninfected control and in vitro M.tb-infected. "
            "Paired design; 2×2 factorial structure with batch effects across collection sites."
        ),
        "sample_groups": (
            "Control (uninfected macrophages) vs M.tb-infected macrophages, "
            "from both healthy donors and active TB patients. "
            "Key comparisons: infected vs uninfected within each donor group; "
            "healthy vs TB patient background effect."
        ),
        "input_file_types": (
            "Raw FASTQ reads deposited in SRA (paired-end, Illumina HiSeq). "
            "Supplementary count matrices provided as tab-separated text files "
            "(gene-level, STAR alignment against GRCh38). "
            "Sample metadata available in GEO series matrix and supplementary Excel sheet."
        ),
        "recommended_start": (
            "Download the supplementary count matrix (GSE148036_raw_counts.txt.gz). "
            "This gives you STAR-generated gene-level counts directly usable in DESeq2/edgeR. "
            "Match sample names to the GEO metadata SOFT file to recover group/batch annotations. "
            "Avoid using the FPKM table for DE analysis."
        ),
        "downstream_workflow": (
            "1. Load count matrix and metadata into R. "
            "2. Filter low-count genes (≥10 counts in ≥50% of samples). "
            "3. Run DESeq2 with design ~ donor + infection_status to account for pairing. "
            "4. LFC shrinkage with apeglm. "
            "5. Pathway enrichment with clusterProfiler (GO, KEGG). "
            "6. GSEA using MSigDB hallmark gene sets to identify broad pathway themes."
        ),
        "reproducibility_notes": (
            "Record the exact STAR and featureCounts version used. "
            "Reference genome and GTF annotation version (Ensembl release 100 vs 104 can "
            "alter gene counts). Donor anonymisation means sample IDs in GEO may differ "
            "from those in the paper — always cross-reference using the Methods section. "
            "This study is relevant to TB immunology and is a good model for paired "
            "macrophage infection designs."
        ),
    },
    {
        "accession": "GSE167232",
        "title": "Single-cell and spatial transcriptomics reveal TB granuloma architecture and cell states",
        "organism": "Homo sapiens / Macaca mulatta",
        "study_question": (
            "What are the spatial organisation and transcriptional cell states within "
            "tuberculosis granulomas, and how do immune and stromal cell niches differ "
            "between contained and progressing granulomas?"
        ),
        "experimental_design": (
            "Combination of bulk RNA-seq (granuloma tissue), scRNA-seq (10x Chromium), "
            "and spatial transcriptomics (Visium) on TB granuloma samples from non-human "
            "primate (rhesus macaque) TB model and human autopsy tissue. "
            "Multiple granuloma types per animal: contained (sterile) vs progressing "
            "(culture-positive). No paired design; batch correction by animal ID."
        ),
        "sample_groups": (
            "Contained granulomas vs progressing granulomas. Within scRNA-seq: "
            "macrophage subsets (foamy, inflammatory, transitional), T cell subsets "
            "(CD4, CD8, regulatory), neutrophils, fibroblasts, epithelioid cells. "
            "Spatial data deconvoluted to infer cell-type composition per spot."
        ),
        "input_file_types": (
            "Bulk RNA-seq: raw counts in supplementary CSV. "
            "scRNA-seq: Cell Ranger output (barcodes, features, matrix) per sample — "
            "can be loaded directly into Seurat or Scanpy. "
            "Spatial Visium: 10x Space Ranger output folders with tissue images and "
            "spot-level count matrices, loadable in Seurat or Squidpy."
        ),
        "recommended_start": (
            "For bulk re-analysis: download supplementary count matrix, apply standard "
            "DESeq2 workflow with ~ granuloma_type as the main effect. "
            "For scRNA-seq: download the processed Seurat RDS object if provided, "
            "or reload raw Cell Ranger outputs and reprocess with Seurat/Scanpy. "
            "Spatial data requires Space Ranger output folder structure — not suitable "
            "for browser-based re-analysis."
        ),
        "downstream_workflow": (
            "Bulk: DESeq2 → pathway enrichment → identify inflammatory vs resolution signatures. "
            "scRNA-seq: QC → normalisation → HVG selection → PCA → UMAP → clustering → "
            "marker gene identification → cell-type annotation → differential abundance. "
            "Spatial: deconvolution with RCTD or cell2location → spatial correlation → "
            "ligand-receptor interaction analysis with CellChat or NicheNet."
        ),
        "reproducibility_notes": (
            "Multi-modal datasets require careful version tracking: Seurat v4 vs v5 "
            "changes the default normalisation workflow. Cell Ranger version affects "
            "barcode calling and alignment. Spatial deconvolution tools are sensitive "
            "to the reference single-cell atlas used. This study illustrates how "
            "complex modern datasets combine multiple data modalities that require "
            "different computational infrastructure — unsuitable for simple web app "
            "re-analysis without HPC resources."
        ),
    },
]

study_options = {f"{cs['accession']} — {cs['title'][:60]}…": i
                 for i, cs in enumerate(CASE_STUDIES)}
selected_label = st.selectbox("Select a study:", list(study_options.keys()))
cs = CASE_STUDIES[study_options[selected_label]]

st.markdown(f"### {cs['accession']}: {cs['title']}")
st.markdown(f"**Organism:** {cs['organism']}")

tabs = st.tabs([
    "🔬 Study Design", "📁 Files & Starting Point",
    "🔧 Workflow", "♻️ Reproducibility Notes"
])

with tabs[0]:
    st.markdown(f"**Research Question:**\n\n{cs['study_question']}")
    st.markdown(f"**Experimental Design:**\n\n{cs['experimental_design']}")
    st.markdown(f"**Sample Groups:**\n\n{cs['sample_groups']}")

with tabs[1]:
    st.markdown(f"**Typical File Types Available:**\n\n{cs['input_file_types']}")
    st.info(f"💡 **Recommended starting file:**\n\n{cs['recommended_start']}", icon="📥")

with tabs[2]:
    st.markdown(f"**Typical Downstream Workflow:**\n\n{cs['downstream_workflow']}")
    st.warning(
        "Public data often need cleaning before use: "
        "sample names in GEO metadata may differ from the paper, "
        "group labels may use abbreviations, and batch/donor annotations "
        "may require manual curation from the supplementary tables.",
        icon="⚠️",
    )

with tabs[3]:
    st.markdown(f"{cs['reproducibility_notes']}")

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 3. COUNTS vs TPM/FPKM EXPLANATION
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("📏 Raw Counts vs TPM / FPKM / CPM")

col_a, col_b = st.columns(2)
with col_a:
    st.markdown("""
**Raw counts** are the number of sequencing reads that align to each gene.
They are integers, directly reflecting sequencing depth and gene expression.

✅ Use raw counts for:
- Differential expression (DESeq2, edgeR, limma-voom)
- Any method that models the count distribution explicitly

❌ Do NOT use raw counts for:
- Cross-sample normalisation comparisons in isolation
  (library size varies between samples)
    """)
with col_b:
    st.markdown("""
**TPM (Transcripts Per Million)**, **FPKM**, and **CPM** are normalised values
that account for library size and (for FPKM/TPM) gene length.

✅ Use normalised values for:
- Visualisation (heatmaps, PCA, boxplots)
- Cross-sample expression comparisons without formal testing
- Single-gene exploratory plots

❌ Do NOT use TPM/FPKM for:
- Input to DESeq2, edgeR — these models require raw counts
- Re-analysis where the original counts are also available
    """)

st.info(
    "**This app expects raw counts.** If you upload TPM or FPKM values, "
    "the statistical models (log-CPM normalisation, t-tests) will still run "
    "but results will be less meaningful, and the data may contain negative "
    "values after log-transformation that trigger validation warnings.",
    icon="⚠️",
)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 4. METADATA HYGIENE
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🧹 Metadata Hygiene")
st.markdown("""
Good metadata is as important as good data. Common issues when re-using
public datasets:

- **Mismatched sample names** — GEO sample IDs (GSM…) rarely match sample
  names in count matrices; always check the supplementary metadata file.
- **Inconsistent group labels** — "control", "ctrl", "Control", "CTRL" are
  four different strings. Standardise before analysis.
- **Missing batch or donor info** — GEO metadata is often incomplete; check
  the Methods and supplementary tables of the paper itself.
- **Ambiguous time points or conditions** — particularly in longitudinal or
  multi-factorial designs; clarify the exact comparison you want before coding groups.
- **Duplicate or swapped samples** — especially in large multi-site studies;
  PCA and sample correlation checks can reveal obvious anomalies.

**Minimum metadata for this app:**
- `groupA` — required; defines comparison groups
- `donor` — optional; enables paired DE analysis
- `batch` — optional; enables batch correction workflow
""")

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 5. REPRODUCIBILITY CHECKLIST (interactive)
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("✅ Reproducibility Checklist")
st.caption(
    "Check off each step as you complete your analysis. "
    "This is a teaching tool — state is not saved between sessions."
)

checklist_items = [
    ("📋 Study context", [
        ("Accession number recorded (e.g. GSE…, E-MTAB-…)", False),
        ("Corresponding paper read and cited", False),
        ("Study design and group definitions noted", False),
        ("Organism and tissue type confirmed", False),
    ]),
    ("📁 Data provenance", [
        ("Raw vs processed files distinguished", False),
        ("Downloaded file checksums verified (if provided)", False),
        ("Gene annotation version noted (e.g. Ensembl release, GTF date)", False),
        ("Reference genome version noted (e.g. GRCh38, hg19)", False),
    ]),
    ("🧬 Metadata", [
        ("Sample metadata inspected and cleaned", False),
        ("Sample names match between count matrix and metadata", False),
        ("Group labels standardised and unambiguous", False),
        ("Batch and donor columns identified if present", False),
    ]),
    ("💻 Computational environment", [
        ("Software versions recorded (R, Python, key packages)", False),
        ("Environment captured (conda env, renv, requirements.txt)", False),
        ("Random seeds set where applicable", False),
        ("Analysis scripts version-controlled (git)", False),
    ]),
    ("📊 Analysis outputs", [
        ("Analysis steps documented with explanatory comments", False),
        ("Output files named consistently and descriptively", False),
        ("Figures traceable to the code that generated them", False),
        ("Filtering and parameter choices justified in notes", False),
        ("Limitations of re-analysis noted", False),
    ]),
]

all_checked = 0
all_total   = 0
for section, items in checklist_items:
    st.markdown(f"**{section}**")
    for label, _ in items:
        checked = st.checkbox(label, key=f"chk_{label}")
        if checked:
            all_checked += 1
        all_total += 1

progress = all_checked / all_total if all_total > 0 else 0
st.progress(progress, text=f"Checklist progress: {all_checked}/{all_total} items")
if all_checked == all_total:
    st.success("🎉 All checklist items complete — your analysis is well-documented!", icon="✅")

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# 6. WEB APP vs LOCAL vs HPC COMPARISON
# ════════════════════════════════════════════════════════════════════════════════

st.subheader("🖥️ Web App vs Local Analysis vs HPC")

comparison_data = {
    "Feature": [
        "Recommended use",
        "Strengths",
        "Limitations",
        "Typical dataset scale",
        "Best fit examples",
    ],
    "This Web App": [
        "Teaching, exploration, small datasets",
        "No setup required; interactive; great for learning core concepts",
        "Limited compute; no file persistence; browser memory constraints; no custom packages",
        "< 5,000 genes × < 100 samples for comfortable performance",
        "Learning RNA-seq, small pilot datasets, teaching demonstrations",
    ],
    "Local Analysis (R/Python)": [
        "Research-grade analysis; full package access",
        "Full DESeq2/edgeR/limma; complete package ecosystem; reproducible scripts",
        "Requires R/Python setup; may be slow on laptop for large datasets",
        "Thousands of genes × hundreds of samples; limited by RAM (typically 8–64 GB)",
        "Standard bulk RNA-seq; scRNA-seq on moderate cohorts; replicating published workflows",
    ],
    "HPC / Server": [
        "Large-scale genomics; alignment pipelines; multi-sample cohorts",
        "Parallel processing; large memory; FASTQ alignment; workflow managers (Snakemake, Nextflow)",
        "Requires HPC access, SLURM knowledge, environment management",
        "Hundreds to thousands of samples; full alignment + quantification pipelines",
        "FASTQ → BAM → counts pipelines; scRNA-seq large cohorts; spatial transcriptomics",
    ],
}

st.dataframe(
    pd.DataFrame(comparison_data).set_index("Feature"),
    use_container_width=True,
)

st.info(
    "This app is designed for teaching and moderate-sized datasets. "
    "For production-grade analysis of real cohorts, use DESeq2/edgeR/limma-voom "
    "locally or on an HPC cluster with proper workflow management.",
    icon="ℹ️",
)

st.divider()

# ════════════════════════════════════════════════════════════════════════════════
# NAVIGATION
# ════════════════════════════════════════════════════════════════════════════════

col_n1, col_n2 = st.columns(2)
with col_n1:
    st.page_link("pages/09_Functional_Enrichment.py",
                 label="← Lesson 9: Functional Enrichment", icon="🧩")
with col_n2:
    st.page_link("pages/11_Single_Cell_RNAseq.py",
                 label="Lesson 11: Single-Cell RNA-seq →", icon="🔬")

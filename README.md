# 🧬 RNA-seq Interactive Training Platform

An interactive, step-by-step Streamlit web application for learning RNA-seq data analysis — from raw counts to biological interpretation. Designed for students, researchers, and bioinformaticians getting started with transcriptomic workflows.

> **Educational disclaimer:** This app is not a production-grade substitute for DESeq2, edgeR, limma-voom, clusterProfiler, fgsea, or full GSEA pipelines. All analyses are implemented as transparent educational models to teach core concepts. For research-grade analysis, use validated R/Python bioinformatics tools.

---

## 🚀 Getting Started

### Run locally

```bash
git clone https://github.com/melodysum/rnaseq-training-platform.git
cd rnaseq-training-platform
pip install -r requirements.txt
streamlit run Home.py
```

### Try online
[rnaseq-training-platform.streamlit.app](https://rnaseq-training-platform.streamlit.app/)

---

## 📚 Course Content

| # | Lesson | Status |
|---|--------|--------|
| 1 | RNA-seq Foundations & Experimental Design | ✅ Available |
| 2 | Quantification, Import & Annotation | ✅ Available |
| 3 | Low-Expression Filtering | ✅ Available |
| 4 | Multiple Testing & FDR | ✅ Available |
| 5 | Exploratory Analysis & PCA | ✅ Available |
| 6 | Batch Correction | ✅ Available |
| 7 | Differential Expression | ✅ Available |
| 8 | Clustering & Heatmaps | ✅ Available |
| 9 | Functional Enrichment (GO / GSEA) | ✅ Available |
| 10 | Public Data & Reproducible Analysis | ✅ Available |
| 11 | Single-Cell RNA-seq (Introductory Module) | ✅ Available |

---

## 📁 Input File Format

The app accepts two CSV files:

### counts.csv — raw gene-level count matrix

- **Genes as rows**, samples as columns
- First column = gene identifiers (used as row index)
- Values must be **raw integer counts** — not TPM, FPKM, or CPM
- At least 2 sample columns and 1 gene row required

```
gene_symbol, D01_control, D01_treatment, D02_control, D02_treatment
IFNA6,       25,          40,            15,          48
CAD,          0,           9,            12,          10
```

### metadata.csv — sample annotation table

- **Samples as rows**, annotations as columns
- Sample names must match the column names in counts.csv exactly (case-sensitive)
- At least 2 samples required

```
sample_name,   groupA,    donor, batch,  sex, age
D01_control,   control,   D01,   batch1, F,   37
D01_treatment, treatment, D01,   batch1, F,   37
D02_control,   control,   D02,   batch2, F,   40
D02_treatment, treatment, D02,   batch2, F,   40
```

---

## 🏷️ Metadata Column Requirements

| Column | Required? | Used for |
|--------|-----------|---------|
| `groupA` | **Required** | Defines comparison groups for all analysis |
| `donor` | Optional | Enables paired DE analysis (Lesson 7) |
| `batch` | Optional | Enables batch correction workflow (Lesson 6) |
| `sex`, `age`, etc. | Optional | Descriptive; displayed in dataset summary |

**Validation rules:**
- Duplicate sample names (in either file) will be rejected
- Duplicate gene identifiers will be rejected
- Non-numeric count values will be rejected
- Negative values will be rejected (counts must be ≥ 0)
- Sample names must match exactly between the two files

---

## 🔬 Feature Scope

### Lessons 1–8 (bulk RNA-seq core)
- Raw count loading and demo data fallback
- Low-expression gene filtering
- Multiple testing / FDR correction demonstration
- PCA-based exploratory analysis
- Simple linear model batch correction
- Educational DE: log-CPM + paired/unpaired t-test + BH FDR
  - Paired DE explicitly aligns samples by donor before testing
  - Falls back to Welch t-test if donor matching is insufficient
- Hierarchical clustering and expression heatmaps

### Lesson 9 — Functional Enrichment
- Built-in educational toy pathway library (10 pathways, real human gene symbols)
- Over-representation analysis (ORA): Fisher's exact test + BH correction
- GSEA-like running-sum enrichment scoring
- Enrichment plot: running sum vs gene rank, pathway hit markers
- Works from session DE results or demo data if no prior DE run

### Lesson 10 — Public Data & Reproducibility
- Real GEO case studies (GSE148036, GSE167232) stored locally
- Interactive reproducibility checklist
- Counts vs TPM/FPKM/CPM explanation
- Web app vs local vs HPC comparison table
- Metadata hygiene guidance

### Lesson 11 — Single-Cell RNA-seq (Introductory)
- Conceptual introduction: cells vs samples, sparsity, UMAP, clustering, marker genes
- Simulated UMAP scatter plot (7 color-coded toy clusters)
- Example marker gene table per cluster
- Bulk vs scRNA-seq structured comparison table
- Explicit "what this app does not implement" section

---

## ⚠️ Educational Limitations

This app is intended for **teaching and moderate-sized datasets**:

- Statistical models are simplified (t-tests on log-CPM, not negative binomial GLMs)
- Not a substitute for DESeq2, edgeR, or limma-voom for differential expression
- Functional enrichment uses toy educational pathways, not real GO/KEGG/MSigDB databases
- Lesson 11 uses simulated UMAP data; no real scRNA-seq pipeline is implemented
- No file persistence between sessions
- Large datasets (>5,000 genes × 200+ samples) may be slow or exceed browser memory
- Runtime depends on deployment resources and dataset size — no specific timings guaranteed
- Interactive clustering and heatmaps should be limited to a reasonable number of selected genes

For production-grade analysis, use:
- **DE:** DESeq2, edgeR, limma-voom
- **Enrichment:** clusterProfiler, fgsea, GSEA desktop
- **scRNA-seq:** Scanpy (Python), Seurat (R)
- **Alignment:** STAR, HISAT2, Salmon + tximeta
- **Workflows:** Snakemake, Nextflow on HPC

---

## 🔬 Why Single-Cell Analysis Is Not Fully Implemented

Lesson 11 is an **introductory conceptual module only**. The following single-cell capabilities are not available in this web app, for the reasons explained below.

| Capability | Why it is not implemented |
|------------|--------------------------|
| Cell-level QC (mitochondrial reads, doublet removal) | Requires per-cell barcode-level data in 10x Cell Ranger format, not CSV |
| Single-cell normalisation (scran, SCTransform) | Specialised methods designed for sparse UMI count distributions; not interchangeable with bulk normalisation |
| Highly variable gene (HVG) selection | Requires per-cell variance modelling across thousands of cells |
| PCA / SNN graph on real scRNA-seq data | Matrices of 20,000 genes × 50,000 cells exceed browser memory limits (~1 GB on Streamlit Cloud) |
| True UMAP on real single-cell matrices | Computationally intensive; real UMAP on large matrices takes minutes and would time out in a web app |
| Leiden / Louvain graph-based clustering | Requires `leidenalg` and `igraph`, which are complex C-extension dependencies incompatible with simple Streamlit deployment |
| Cluster-level differential expression | Depends on completed clustering; not meaningful without real cell-type assignments |
| Cell-type annotation from marker genes | Requires a reference atlas or curated marker database not bundled with this app |
| Multi-sample integration (Harmony, scVI, Seurat CCA) | Requires GPU or large-memory server; scVI depends on PyTorch |
| Trajectory / pseudotime analysis (Monocle, scVelo) | Requires RNA velocity or pseudotime-specific data structures beyond count matrices |
| AnnData (.h5ad) / Seurat (.rds) file support | Binary HDF5 formats; files are typically several GB and cannot be uploaded via browser |
| 10x Cell Ranger output parsing | Requires reading sparse matrix market format (barcodes / features / matrix triplet), not CSV |

**In short:** scRNA-seq analysis requires specialised file formats, large memory (8–64 GB RAM), complex library dependencies, and sometimes GPU resources. These constraints make it unsuitable for a lightweight browser-based teaching app.

**For real single-cell analysis, use:**
- [Scanpy](https://scanpy.readthedocs.io/) (Python) — standard single-cell toolkit
- [Seurat](https://satijalab.org/seurat/) (R) — widely used in the field
- [Bioconductor scran / scater](https://bioconductor.org/books/release/OSCA/) (R) — rigorous statistical framework
- Run on a local machine with ≥16 GB RAM, or on an HPC cluster

---

## 🗺️ Roadmap / Future Extensions

- Lessons 12+: trajectory analysis concepts, multi-omics integration overview
- Real GO/KEGG enrichment via offline MSigDB GMT file support
- DESeq2-style negative binomial modelling (via rpy2 or pydeseq2)
- Better support for multi-group / multi-factor designs
- Improved spatial transcriptomics conceptual module

---

## 📂 Project Structure

```
rnaseq-training-platform/
├── Home.py                          # Landing page + data upload
├── requirements.txt
├── data/
│   ├── counts.csv                   # Built-in demo count matrix (10,000 genes × 40 samples)
│   └── metadata.csv                 # Built-in demo metadata (paired design, 20 donors)
├── pages/
│   ├── 01_RNAseq_Foundations.py
│   ├── 02_Quantification_Import_Annotation.py
│   ├── 03_Filtering.py
│   ├── 04_FDR.py
│   ├── 05_Exploratory_Analysis_PCA.py
│   ├── 06_Batch_Correction.py
│   ├── 07_Differential_Expression.py
│   ├── 08_Clustering_Heatmaps.py
│   ├── 09_Functional_Enrichment.py  # ORA + GSEA-like enrichment
│   ├── 10_Public_Data.py            # GEO case studies + reproducibility
│   └── 11_Single_Cell_RNAseq.py    # scRNA-seq conceptual intro
└── utils/
    ├── data_loader.py               # Loading, validation, session init
    ├── de_analysis.py               # Educational DE workflow
    ├── enrichment_utils.py          # Toy gene sets, ORA, GSEA-like
    ├── batch_effects.py
    ├── clustering_utils.py
    ├── exploration.py
    ├── fdr_demo.py
    ├── filtering.py
    ├── pca_utils.py
    └── simulation.py
```

---

## 📄 License

MIT License. See LICENSE file.

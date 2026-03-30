# RNA-seq Interactive Training Platform

[![Live app](https://img.shields.io/badge/Live%20App-Streamlit-ff4b4b?logo=streamlit&logoColor=white)](https://rnaseq-training-platform.streamlit.app/)
[![GitHub repo](https://img.shields.io/badge/GitHub-melodysum%2Frnaseq--training--platform-181717?logo=github)](https://github.com/melodysum/rnaseq-training-platform)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](#installation)
[![Status](https://img.shields.io/badge/Status-Teaching%20Prototype-6f42c1)](#important-note)
[![Curriculum](https://img.shields.io/badge/Curriculum-Updated%20to%20Lesson%208-0a7ea4)](#course-roadmap)

An interactive **Streamlit teaching platform** for learning core RNA-seq analysis concepts through guided visual exploration.

This project is designed as a **student-friendly educational website**, not a production-grade RNA-seq analysis pipeline. It helps learners move from conceptual foundations to statistical interpretation through modular, interactive lessons.

## Live app

**Streamlit app:** https://rnaseq-training-platform.streamlit.app/  
**Source code:** https://github.com/melodysum/rnaseq-training-platform

---

## What this project is

The platform turns static RNA-seq teaching material into an interactive curriculum with:

- concept-first lessons
- guided visual explanations
- uploadable user datasets
- built-in demo data fallback
- progressive navigation from foundations to downstream interpretation

It is intended for:

- students learning bulk RNA-seq analysis
- wet-lab researchers who want clearer intuition for transcriptomics workflows
- teaching demonstrations and classroom-style guided exploration

---

## Current lesson coverage

The app is currently updated through **Lesson 8**.

### Implemented lessons

1. **RNA-seq Foundations & Experimental Design**  
   Sequencing basics, FASTQ structure, replication, depth, and study design logic.

2. **Quantification, Import & Annotation**  
   From FASTQ to quantified expression, transcript-to-gene aggregation, metadata, and annotation.

3. **Low-Expression Filtering**  
   Why low-count genes are filtered and how threshold choices affect retention.

4. **Multiple Testing & FDR**  
   BH correction, testing burden, and the relationship between filtering and adjusted significance.

5. **Exploratory Analysis & PCA**  
   Sample structure, variance, metadata-aware PCA, and outlier exploration.

6. **Batch Correction**  
   Technical variation, confounding, correction logic, and downstream consequences.

7. **Differential Expression**  
   Group comparison, log fold change, adjusted p-values, and result interpretation.

8. **Clustering & Heatmaps**  
   Gene/sample clustering, row scaling, heatmap logic, and gene module exploration.

### Planned lessons

9. **Functional Enrichment (GO / GSEA)**  
   Pathway-level interpretation of DE results.

10. **Public Data & Reproducible Analysis**  
   Public transcriptomics resources, HDF5-style data access, and reproducible workflows.

11. **Introduction to Single-cell RNA-seq**  
   10x, UMI, droplets, scRNA-seq QC, clustering, UMAP, and cell type annotation.

---

## Course roadmap

| Lesson | Topic | Status | Main focus |
|---|---|---:|---|
| 1 | RNA-seq Foundations & Experimental Design | ✅ Implemented | Sequencing concepts, replicates, design trade-offs |
| 2 | Quantification, Import & Annotation | ✅ Implemented | Quantification logic, TPM, annotation, metadata |
| 3 | Low-Expression Filtering | ✅ Implemented | Thresholding, retained genes, filtering rationale |
| 4 | Multiple Testing & FDR | ✅ Implemented | BH correction, adjusted p-values, testing burden |
| 5 | Exploratory Analysis & PCA | ✅ Implemented | Sample structure, variance, outliers |
| 6 | Batch Correction | ✅ Implemented | Batch effects, confounding, correction |
| 7 | Differential Expression | ✅ Implemented | Group comparison, volcano logic, DE outputs |
| 8 | Clustering & Heatmaps | ✅ Implemented | Heatmaps, clustering, modules |
| 9 | Functional Enrichment (GO / GSEA) | 🟡 Planned | ORA, GSEA, pathway interpretation |
| 10 | Public Data & Reproducible Analysis | 🟡 Planned | Public resources, HDF5, reproducibility |
| 11 | Introduction to Single-cell RNA-seq | 🟡 Planned | 10x, UMI, QC, UMAP, annotation |

---

## Screenshots

> To make these render on GitHub, place screenshots in `docs/screenshots/` inside the repository.

### Home / course overview

![Home / course overview](docs/screenshots/home-course-overview.png)

### Suggested additional screenshots

You can add more images later and reference them here, for example:

```md
![Lesson 3 — Filtering](docs/screenshots/lesson3-filtering.png)
![Lesson 6 — Batch Correction](docs/screenshots/lesson6-batch-correction.png)
![Lesson 8 — Clustering & Heatmaps](docs/screenshots/lesson8-clustering-heatmaps.png)
```

---

## Important note

This app is built for **teaching and exploration**.

Several modules include **simplified educational implementations** so students can understand the logic of RNA-seq analysis without requiring a full production bioinformatics environment.

It is **not a substitute** for publication-grade workflows such as:

- **DESeq2**
- **edgeR**
- **limma-voom**
- specialised enrichment frameworks
- full single-cell pipelines such as Seurat / Scanpy end-to-end analyses

Use this platform to understand concepts, inspect patterns, and build intuition — not to generate final analytical results for publication.

---

## Input file format

The app accepts two CSV files.

### `counts.csv`
Raw count matrix.

Requirements:
- first column = gene identifier column
- remaining columns = sample names
- values = raw counts

Example:

```csv
gene_symbol,D01_control,D01_treatment,D02_control,D02_treatment
IFNA6,25,40,15,48
CAD,0,9,12,10
```

### `metadata.csv`
One row per sample.

Requirements:
- sample names must match the columns in `counts.csv`
- must contain metadata suitable for downstream grouping / interpretation

Recommended columns:
- `sample_name`
- `groupA`
- `donor`
- `batch`
- `sex`
- `age`

Example:

```csv
sample_name,groupA,donor,batch,sex,age
D01_control,control,D01,batch1,F,37
D01_treatment,treatment,D01,batch1,F,37
```

### Validation behaviour

The app checks that:
- uploaded samples match between counts and metadata
- required fields are present
- lesson modules can reuse the uploaded dataset consistently through session state

If no files are uploaded, the app falls back to the built-in demo dataset in `data/`.

---

## Project structure

```text
rnaseq-training-platform/
├── Home.py
├── requirements.txt
├── data/
│   ├── counts.csv
│   └── metadata.csv
├── pages/
│   ├── 01_RNAseq_Foundations.py
│   ├── 02_Quantification_Import_Annotation.py
│   ├── 03_Filtering.py
│   ├── 04_FDR.py
│   ├── 05_Exploratory_Analysis_PCA.py
│   ├── 06_Batch_Correction.py
│   ├── 07_Differential_Expression.py
│   ├── 08_Clustering_Heatmaps.py
│   ├── 09_Functional_Enrichment.py
│   └── 10_Public_Data.py
└── utils/
    ├── data_loader.py
    ├── filtering.py
    ├── fdr_demo.py
    ├── pca_utils.py
    ├── de_analysis.py
    ├── batch_effects.py
    ├── simulation.py
    ├── exploration.py
    └── clustering_utils.py
```

### File overview

- `Home.py`  
  Course landing page, lesson navigation, data upload, and dataset summary.

- `pages/01_RNAseq_Foundations.py`  
  Foundational lesson on sequencing concepts and experimental design.

- `pages/02_Quantification_Import_Annotation.py`  
  Educational lesson on quantification, annotation, and metadata structure.

- `pages/03_Filtering.py`  
  Interactive lesson for low-expression filtering.

- `pages/04_FDR.py`  
  Educational lesson for BH/FDR and downstream significance behaviour.

- `pages/05_Exploratory_Analysis_PCA.py`  
  PCA, variance structure, and outlier exploration.

- `pages/06_Batch_Correction.py`  
  Batch-effect detection and correction demonstrations.

- `pages/07_Differential_Expression.py`  
  Teaching-oriented DE analysis and result inspection.

- `pages/08_Clustering_Heatmaps.py`  
  Gene/sample clustering, heatmaps, and module exploration.

- `pages/09_Functional_Enrichment.py`  
  Planned pathway interpretation module.

- `pages/10_Public_Data.py`  
  Planned public-data and reproducibility module.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/melodysum/rnaseq-training-platform.git
cd rnaseq-training-platform
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate   # macOS / Linux
```

On Windows:

```bash
.venv\Scripts\activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the app

```bash
streamlit run Home.py
```

---

## Dependencies

Main packages used in the current app:

- `streamlit`
- `pandas`
- `numpy`
- `scipy`
- `statsmodels`
- `plotly`

---

## How to use

### Option 1 — Use demo data
Launch the app and keep the built-in dataset selected.

### Option 2 — Upload your own data
On the Home page:
1. upload `counts.csv`
2. upload `metadata.csv`
3. move through lessons in sequence or jump directly to a module of interest

### Suggested learning flow
1. Start with **Lesson 1–2** for conceptual foundations
2. Move to **Filtering** and **FDR**
3. Use **PCA** and **Batch Correction** to inspect sample structure
4. Run **Differential Expression**
5. Explore **Clustering & Heatmaps**
6. Follow future updates for enrichment, public data, and single-cell modules

---

## Current limitations

- This is still a teaching prototype rather than a production analysis suite
- Several analytical modules use simplified educational logic for speed and clarity
- Lessons 9–11 are planned but not yet fully implemented in the current version
- Input validation is tailored to the current expected CSV structure
- Large interactive matrices may still require careful truncation for performance

---

## Planned extensions

Planned future development includes:

- Functional enrichment (GO / GSEA)
- Public RNA-seq data reuse and reproducibility modules
- Introductory single-cell RNA-seq teaching module
- richer visual teaching components
- improved screenshot coverage and lesson-level documentation

---

## Use case

This platform is especially useful for:

- teaching RNA-seq concepts step by step
- helping students understand how analysis choices change results
- demonstrating the relationship between preprocessing, statistics, and interpretation
- creating guided transcriptomics learning experiences in a web interface

---

## Acknowledgement

This prototype was developed as part of a broader effort to transform RNA-seq teaching material into a more dynamic, modular, and interactive training website.

---

## License

Add a license here if you plan to make the project openly reusable.

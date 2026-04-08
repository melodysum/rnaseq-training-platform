"""
utils/enrichment_utils.py
--------------------------
Educational enrichment analysis utilities for the RNA-seq Training Platform.

Provides:
  - TOY_GENE_SETS  : built-in educational pathway library (real human gene symbols)
  - run_ora()      : over-representation analysis (Fisher's exact + BH correction)
  - run_gsea_like(): running-sum enrichment scoring (KS/GSEA-style)
  - demo_ranked_genes(): generate a ranked gene list from built-in demo data

NOT a substitute for clusterProfiler, fgsea, GSEA, or MSigDB-based analysis.
All pathways are educational toy sets for teaching purposes only.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ═══════════════════════════════════════════════════════════════════════════════
# TOY GENE SET LIBRARY
# Real human gene symbols; loosely modelled on GO / MSigDB-style categories.
# Each set contains 15–25 genes chosen to be biologically recognisable.
# ═══════════════════════════════════════════════════════════════════════════════

TOY_GENE_SETS = {

    # ── 1. Interferon response ────────────────────────────────────────────────
    # Anchored to canonical ISG and IRF/STAT signalling genes
    "INTERFERON_RESPONSE": [
        "IFNA1", "IFNB1", "IFNG", "IRF1", "IRF3", "IRF7", "IRF9",
        "STAT1", "STAT2", "MX1", "MX2", "OAS1", "OAS2", "ISG15", "ISG20",
        "IFIT1", "IFIT2", "IFIT3", "RSAD2", "CXCL10", "DDX58", "IFITM1",
        "IFITM3", "BST2", "HERC5",
    ],

    # ── 2. Inflammation / TNF / NF-kB ─────────────────────────────────────────
    # Core inflammatory cytokines and NF-kB pathway members
    "INFLAMMATION_NFKB": [
        "TNF", "IL1A", "IL1B", "IL6", "CXCL8", "CXCL1", "CXCL2",
        "CCL2", "CCL5", "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB",
        "PTGS2", "MMP9", "VCAM1", "ICAM1", "TNFAIP3", "CXCL3",
        "IL18", "NLRP3", "CASP1",
    ],

    # ── 3. Cell cycle / mitosis ────────────────────────────────────────────────
    # Anchored to known mitotic and S-phase regulators
    "CELL_CYCLE_MITOSIS": [
        "MKI67", "PCNA", "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDC25C",
        "CDK1", "AURKA", "AURKB", "BUB1", "BUB1B", "PLK1", "TOP2A",
        "UBE2C", "CENPE", "CENPF", "KIF11", "KIF20A", "NUSAP1",
        "CDKN2A", "E2F1", "MCM2",
    ],

    # ── 4. Mitochondrial / oxidative phosphorylation ──────────────────────────
    # Complex I–V subunits and mitochondrial biogenesis factors
    "OXPHOS_MITOCHONDRIA": [
        "NDUFA1", "NDUFA2", "NDUFS1", "NDUFB1", "SDHA", "SDHB",
        "UQCRC1", "UQCRB", "COX4I1", "COX5A", "ATP5A1", "ATP5B",
        "CYCS", "TFAM", "VDAC1", "TOMM20", "TIMM23", "SLC25A4",
        "IDH2", "FH", "MDH2", "SUCLA2",
    ],

    # ── 5. Ribosome / translation ──────────────────────────────────────────────
    # Ribosomal proteins and translation initiation factors
    "RIBOSOME_TRANSLATION": [
        "RPS3", "RPS4X", "RPS6", "RPS7", "RPS11", "RPS14",
        "RPL3", "RPL4", "RPL6", "RPL10", "RPL13", "RPL18",
        "RPL23", "RPL26", "EIF4E", "EIF4G1", "EIF2S1", "PABPC1",
        "EIF3A", "EIF1AX", "RPLP0",
    ],

    # ── 6. Epithelial markers ──────────────────────────────────────────────────
    # Tight junction, keratin, and epithelial identity genes
    "EPITHELIAL_MARKERS": [
        "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7",
        "OCLN", "CLDN3", "CLDN4", "TJP1", "MUC1", "ESRP1",
        "ESRP2", "RAB25", "GRHL2", "OVOL2", "ST14", "DSP",
        "PKP3", "CLDN7",
    ],

    # ── 7. Apoptosis / stress response ────────────────────────────────────────
    # Intrinsic/extrinsic apoptosis and DNA damage response
    "APOPTOSIS_STRESS": [
        "TP53", "BAX", "BAD", "BCL2", "BCL2L1", "CASP3", "CASP7",
        "CASP8", "CASP9", "PARP1", "CYCS", "APAF1", "BID",
        "TNFRSF10A", "TNFRSF10B", "GADD45A", "GADD45B", "CDKN1A",
        "MDM2", "DDIT3", "ATM", "CHEK2",
    ],

    # ── 8. Antigen presentation / immune activation ───────────────────────────
    # MHC class I & II, antigen processing, and immune co-stimulation
    "ANTIGEN_PRESENTATION": [
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
        "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "B2M",
        "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "CD74",
        "CIITA", "NLRC5", "CD80", "CD86",
    ],

    # ── 9. Hypoxia response ────────────────────────────────────────────────────
    # HIF targets and glycolytic shift genes
    "HYPOXIA_RESPONSE": [
        "HIF1A", "EPAS1", "VEGFA", "VEGFC", "SLC2A1", "LDHA",
        "PGAM1", "ENO1", "ALDOA", "PKM", "TPI1", "BNIP3",
        "BNIP3L", "HMOX1", "CA9", "EGLN1", "EGLN3", "PDK1",
        "PFKL", "PGK1",
    ],

    # ── 10. Lipid / metabolic reprogramming ───────────────────────────────────
    # Fatty acid synthesis, cholesterol, and nuclear receptor targets
    "LIPID_METABOLISM": [
        "FASN", "ACACA", "ACACB", "SCD", "ELOVL1", "ELOVL6",
        "DGAT1", "DGAT2", "PPARA", "PPARG", "PPARD", "SREBF1",
        "SREBF2", "INSIG1", "INSIG2", "HMGCR", "LDLR", "ACLY",
        "ACSL4", "FADS1",
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DEMO DATA HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def demo_ranked_genes(counts: pd.DataFrame, metadata: pd.DataFrame,
                      group_col: str = "groupA") -> pd.DataFrame:
    """
    Generate a ranked gene list from built-in demo data using log2FC
    between the first two groups found in metadata[group_col].

    Returns a DataFrame with columns: gene, log2FC, mean_g1, mean_g2
    sorted by log2FC descending (for GSEA-like input).
    """
    groups = metadata[group_col].dropna().unique()
    if len(groups) < 2:
        raise ValueError("Need at least 2 groups in metadata for demo ranking.")

    g1, g2 = sorted(groups)[:2]
    s1 = metadata[metadata[group_col] == g1].index.tolist()
    s2 = metadata[metadata[group_col] == g2].index.tolist()
    s1 = [s for s in s1 if s in counts.columns]
    s2 = [s for s in s2 if s in counts.columns]

    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    lcpm = np.log2(cpm + 1)

    m1 = lcpm[s1].mean(axis=1)
    m2 = lcpm[s2].mean(axis=1)
    lfc = m2 - m1

    df = pd.DataFrame({
        "gene":    counts.index,
        "log2FC":  lfc.values,
        "mean_g1": m1.values,
        "mean_g2": m2.values,
    })
    df["group1"] = g1
    df["group2"] = g2
    return df.sort_values("log2FC", ascending=False).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
# OVER-REPRESENTATION ANALYSIS (ORA)
# ═══════════════════════════════════════════════════════════════════════════════

def run_ora(gene_list: list,
            gene_sets: dict = None,
            universe: list = None) -> pd.DataFrame:
    """
    Over-representation analysis using Fisher's exact test + BH correction.

    Parameters
    ----------
    gene_list : list of gene symbols (query set)
    gene_sets : dict {pathway_name: [gene_symbols]}; defaults to TOY_GENE_SETS
    universe  : background gene list; defaults to all unique genes in gene_sets

    Returns
    -------
    DataFrame with columns:
        pathway, pathway_size, query_size, overlap, gene_ratio,
        pvalue, padj, overlap_genes
    sorted by padj ascending.
    """
    if gene_sets is None:
        gene_sets = TOY_GENE_SETS

    # Build universe from all pathway genes if not supplied
    if universe is None or len(universe) == 0:
        universe = list({g for gs in gene_sets.values() for g in gs})

    universe_set = set(universe)
    query_set    = set(g for g in gene_list if g in universe_set)
    N = len(universe_set)   # total background genes
    K = len(query_set)      # query genes in background

    rows = []
    for pathway, members in gene_sets.items():
        pathway_set   = set(members) & universe_set
        M             = len(pathway_set)           # pathway size in background
        overlap_genes = sorted(query_set & pathway_set)
        k             = len(overlap_genes)         # overlap

        # Fisher's exact test (one-sided: enrichment)
        # Contingency table:
        #             in_pathway | not_in_pathway
        # in_query:      k       |   K - k
        # not_query:   M - k     |   N - M - (K - k)
        table = [
            [k,         K - k],
            [M - k,     N - M - (K - k)],
        ]
        _, pval = stats.fisher_exact(table, alternative="greater")

        rows.append({
            "pathway":      pathway,
            "pathway_size": M,
            "query_size":   K,
            "overlap":      k,
            "gene_ratio":   round(k / M, 4) if M > 0 else 0,
            "pvalue":       pval,
            "overlap_genes": ", ".join(overlap_genes),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["padj"] = padj
    df = df.sort_values("padj").reset_index(drop=True)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# GSEA-LIKE RUNNING SUM ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════

def run_gsea_like(ranked_genes: pd.Series,
                  gene_sets: dict = None,
                  lfc_weights: pd.Series = None) -> pd.DataFrame:
    """
    Running-sum enrichment scoring, conceptually similar to KS/GSEA.

    Parameters
    ----------
    ranked_genes : pd.Series of gene names in ranked order (best to worst).
                   Index should be integer rank position.
    gene_sets    : dict {pathway_name: [gene_symbols]}; defaults to TOY_GENE_SETS
    lfc_weights  : pd.Series {gene_name: abs(log2FC)} for weighted stepping.
                   If None, uses unweighted ±1/N stepping.

    Returns
    -------
    DataFrame with columns:
        pathway, ES, direction, n_hits, pathway_size
    sorted by abs(ES) descending.
    """
    if gene_sets is None:
        gene_sets = TOY_GENE_SETS

    gene_list = list(ranked_genes)
    N         = len(gene_list)
    results   = []

    for pathway, members in gene_sets.items():
        pathway_set = set(members)
        hits        = [g in pathway_set for g in gene_list]
        n_hits      = sum(hits)

        if n_hits == 0:
            results.append({
                "pathway":      pathway,
                "ES":           0.0,
                "direction":    "—",
                "n_hits":       0,
                "pathway_size": len(members),
            })
            continue

        n_miss = N - n_hits

        # Weight for hit positions: abs(log2FC) / sum_abs_lfc_hits
        if lfc_weights is not None:
            lfc_in_path = [abs(lfc_weights.get(g, 0.0)) for g, h in zip(gene_list, hits) if h]
            total_weight = sum(lfc_in_path) or 1.0
            hit_step  = [abs(lfc_weights.get(g, 0.0)) / total_weight if h else 0
                         for g, h in zip(gene_list, hits)]
        else:
            hit_step = [1.0 / n_hits if h else 0 for h in hits]

        miss_step = 1.0 / n_miss if n_miss > 0 else 0.0

        # Walk down the ranked list accumulating the running sum
        running_sum = []
        es = 0.0
        for i, is_hit in enumerate(hits):
            if is_hit:
                es += hit_step[i]
            else:
                es -= miss_step
            running_sum.append(es)

        # Enrichment score = maximum deviation from zero
        max_pos = max(running_sum)
        max_neg = min(running_sum)
        if abs(max_pos) >= abs(max_neg):
            enrichment_score = max_pos
            direction = "Positive (up)"
        else:
            enrichment_score = max_neg
            direction = "Negative (down)"

        results.append({
            "pathway":      pathway,
            "ES":           round(enrichment_score, 4),
            "direction":    direction,
            "n_hits":       n_hits,
            "pathway_size": len(members),
            "_running_sum": running_sum,   # kept for plotting, dropped in return
        })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    df = df.sort_values("ES", key=abs, ascending=False).reset_index(drop=True)
    return df


def get_running_sum(ranked_genes: pd.Series,
                    pathway_genes: list,
                    lfc_weights: pd.Series = None) -> list:
    """
    Return the running-sum vector for a single pathway (for enrichment plots).

    Parameters
    ----------
    ranked_genes   : ordered list/Series of gene names
    pathway_genes  : genes belonging to this pathway
    lfc_weights    : abs(log2FC) per gene for weighted stepping

    Returns
    -------
    list of float — running enrichment score at each rank position
    """
    gene_list   = list(ranked_genes)
    N           = len(gene_list)
    pathway_set = set(pathway_genes)
    hits        = [g in pathway_set for g in gene_list]
    n_hits      = sum(hits)

    if n_hits == 0:
        return [0.0] * N

    n_miss = N - n_hits

    if lfc_weights is not None:
        lfc_in_path  = [abs(lfc_weights.get(g, 0.0)) for g, h in zip(gene_list, hits) if h]
        total_weight = sum(lfc_in_path) or 1.0
        hit_step     = [abs(lfc_weights.get(g, 0.0)) / total_weight if h else 0
                        for g, h in zip(gene_list, hits)]
    else:
        hit_step = [1.0 / n_hits if h else 0 for h in hits]

    miss_step   = 1.0 / n_miss if n_miss > 0 else 0.0
    running_sum = []
    es = 0.0
    for i, is_hit in enumerate(hits):
        es += hit_step[i] if is_hit else -miss_step
        running_sum.append(es)

    return running_sum

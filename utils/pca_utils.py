"""
utils/pca_utils.py
------------------
PCA helper functions for the Batch Correction lesson.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(expr_matrix: pd.DataFrame,
            n_components: int = 10,
            scale: bool = True,
            top_var_genes: int = None) -> tuple:
    """
    Run PCA on an expression matrix (genes × samples).

    Parameters
    ----------
    expr_matrix   : genes × samples DataFrame
    n_components  : number of PCs to compute
    scale         : whether to scale features (samples) before PCA
    top_var_genes : if set, subset to top N most variable genes first

    Returns
    -------
    scores     : pd.DataFrame, shape (samples, n_components)
    explained  : np.array of variance explained ratios
    loadings   : pd.DataFrame, shape (genes, n_components)
    """
    mat = expr_matrix.copy()

    if top_var_genes and top_var_genes < len(mat):
        gene_var = mat.var(axis=1)
        top_idx  = gene_var.nlargest(top_var_genes).index
        mat      = mat.loc[top_idx]

    # PCA is performed on samples, so transpose: rows = samples
    X = mat.T.values
    n_comp = min(n_components, X.shape[0], X.shape[1])

    if scale:
        X = StandardScaler().fit_transform(X)

    pca = PCA(n_components=n_comp)
    scores_arr = pca.fit_transform(X)

    pc_names = [f"PC{i+1}" for i in range(n_comp)]
    scores   = pd.DataFrame(scores_arr, index=mat.columns, columns=pc_names)
    loadings = pd.DataFrame(
        pca.components_.T,
        index=mat.index,
        columns=pc_names,
    )
    return scores, pca.explained_variance_ratio_, loadings


def pca_plot_df(scores: pd.DataFrame,
                metadata: pd.DataFrame,
                pc_x: str = "PC1",
                pc_y: str = "PC2") -> pd.DataFrame:
    """
    Merge PCA scores with metadata for plotting.
    Returns a DataFrame ready for px.scatter.
    """
    df = scores[[pc_x, pc_y]].copy()
    df.index.name = "sample"
    df = df.reset_index()
    if metadata is not None:
        meta = metadata.reset_index().rename(columns={metadata.index.name or "index": "sample"})
        df = df.merge(meta, on="sample", how="left")
    return df

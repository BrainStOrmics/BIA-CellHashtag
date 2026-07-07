"""Find differentially expressed genes per cluster."""
from typing import Any


def find_markers(adata: Any, groupby: str, method: str = "wilcoxon", top_n: int = 20) -> dict[str, list[str]]:
    import scanpy as sc
    import numpy as np

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not in adata.obs")

    if "rank_genes_groups" not in adata.uns or adata.uns.get("_rank_genes_groups_key") != groupby:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)

    clusters = adata.obs[groupby].astype(str).unique()
    markers = {}
    for cluster in clusters:
        scores = adata.uns["rank_genes_groups"]["scores"][cluster]
        names = adata.uns["rank_genes_groups"]["names"][cluster]
        order = np.argsort(scores)[::-1]
        markers[str(cluster)] = [names[i] for i in order[:top_n]]

    return markers

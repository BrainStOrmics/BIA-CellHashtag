"""Find consensus markers across multiple DEG methods."""
from typing import Any


def consensus_markers(
    adata: Any,
    groupby: str,
    methods: tuple[str, ...] = ("wilcoxon", "t-test"),
    min_methods: int = 2,
    top_n: int = 20,
) -> dict[str, list[str]]:
    import scanpy as sc
    import numpy as np
    from collections import Counter

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not in adata.obs")

    clusters = adata.obs[groupby].astype(str).unique()
    all_markers: dict[str, list[set[str]]] = {}

    for method in methods:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
        for cluster in clusters:
            scores = adata.uns["rank_genes_groups"]["scores"][cluster]
            names = adata.uns["rank_genes_groups"]["names"][cluster]
            order = np.argsort(scores)[::-1]
            cluster_id = str(cluster)
            if cluster_id not in all_markers:
                all_markers[cluster_id] = []
            all_markers[cluster_id].append(set(names[i] for i in order[:top_n]))

    consensus = {}
    for cluster_id, marker_sets in all_markers.items():
        if len(marker_sets) < min_methods:
            consensus[cluster_id] = []
            continue
        counts = Counter()
        for s in marker_sets:
            counts.update(s)
        consensus[cluster_id] = [g for g, c in counts.most_common(top_n) if c >= min_methods]

    return consensus

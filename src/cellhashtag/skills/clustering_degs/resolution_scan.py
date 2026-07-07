"""Scan clustering resolution and return quality scores."""
from typing import Any


def resolution_scan(
    adata: Any,
    resolutions: list[float],
    cluster_key: str = "leiden",
    quality_fn: str = "silhouette",
) -> dict[float, float]:
    import scanpy as sc
    from cellhashtag.skills.clustering_quality import silhouette_check, modularity_check

    qfns = {"silhouette": lambda a, k: silhouette_check(a, k)["asw"],
            "modularity": lambda a, k: modularity_check(a, k)["modularity"]}
    if quality_fn not in qfns:
        raise ValueError(f"quality_fn must be one of {list(qfns.keys())}")

    scores = {}
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=cluster_key)
        scores[res] = qfns[quality_fn](adata, cluster_key)

    return scores

"""
Silhouette 轮廓系数检查。

使用 sklearn 的 silhouette_score 评估聚类分离度。
ASW (Average Silhouette Width) ∈ [-1, 1]:
  > 0.5  → 聚类结构清晰
  0.25-0.5 → 中等
  < 0.25 → 聚类可能不清晰
"""

from typing import Any

try:
    import numpy as np
    from sklearn.metrics import silhouette_score, silhouette_samples
except ImportError:
    np = None
    silhouette_score = None
    silhouette_samples = None


def silhouette_check(
    adata: Any,
    cluster_key: str = "leiden",
    use_rep: str = "X_pca",
) -> dict:
    """
    计算 Average Silhouette Width (ASW) 和 per-cluster silhouette。

    Args:
        adata: AnnData 对象。
        cluster_key: adata.obs 中的聚类列名。
        use_rep: 用于计算距离的 embedding 键名。

    Returns:
        {
            "tool": "silhouette",
            "asw": float,              # 全局 ASW
            "per_cluster": dict,        # {cluster_id: asw}
            "recommendation": str,      # "ok" | "adjust" | "maybe_overcluster"
            "details": str,             # 人类可读摘要
        }
    """
    if silhouette_score is None:
        raise ImportError(
            "sklearn is required. Install with: pip install scikit-learn"
        )

    if use_rep not in adata.obsm:
        raise ValueError(f"Embedding '{use_rep}' not found in adata.obsm")

    labels = adata.obs[cluster_key].astype(str)
    embedding = adata.obsm[use_rep]

    # 全局 ASW
    asw = float(silhouette_score(embedding, labels, metric="euclidean"))

    # Per-cluster silhouette
    per_cell_scores = silhouette_samples(embedding, labels, metric="euclidean")
    per_cluster = {}
    for cluster in adata.obs[cluster_key].unique():
        mask = adata.obs[cluster_key] == cluster
        per_cluster[str(cluster)] = float(per_cell_scores[mask].mean())

    # 决策
    n_clusters = labels.nunique()
    if asw < 0.25:
        recommendation = "adjust"
    elif asw > 0.7 and n_clusters > 20:
        recommendation = "maybe_overcluster"
    else:
        recommendation = "ok"

    return {
        "tool": "silhouette",
        "asw": round(asw, 4),
        "per_cluster": {k: round(v, 4) for k, v in per_cluster.items()},
        "recommendation": recommendation,
        "details": f"ASW={asw:.3f}, {n_clusters} clusters",
    }

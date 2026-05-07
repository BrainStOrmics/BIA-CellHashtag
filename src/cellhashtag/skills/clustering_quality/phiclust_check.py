"""
Phiclust 检查 — 基于随机矩阵理论检测 cluster 内隐藏亚群。

Phiclust 是 R 包，通过 rpy2 调用。
φ_clust 高 → cluster 内还有有意义的生物学变异 → 需要提高 resolution
φ_clust 低 → cluster 内部没有结构 → resolution 合适

论文: Mircea et al., Genome Biology (2022)
GitHub: https://github.com/semraulab/phiclust
"""

from typing import Any


def phiclust_check(
    adata: Any,
    cluster_key: str = "leiden",
    threshold: float = 0.5,
) -> dict:
    """
    对每个 cluster 计算 Phiclust 分数。

    Args:
        adata: AnnData 对象。
        cluster_key: adata.obs 中的聚类列名。
        threshold: φ 分数阈值，高于此值认为有隐藏亚群。

    Returns:
        {
            "tool": "phiclust",
            "per_cluster_scores": dict,   # {cluster_id: phi_score}
            "has_substructure": dict,      # {cluster_id: bool}
            "recommendation": str,         # "increase_resolution" | "ok" | "decrease_resolution"
            "details": str,
        }
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
    except ImportError:
        # rpy2 不可用时返回降级结果
        return {
            "tool": "phiclust",
            "per_cluster_scores": {},
            "has_substructure": {},
            "recommendation": "ok",
            "details": "Phiclust unavailable (rpy2 not installed). Skipping R-based check.",
            "warning": "Install rpy2 and the phiclust R package to enable this check.",
        }

    numpy2ri.activate()
    pandas2ri.activate()

    try:
        # 检查 phiclust R 包是否安装
        ro.r('if (!requireNamespace("phiclust", quietly = TRUE)) stop("phiclust R package not installed")')
    except Exception:
        return {
            "tool": "phiclust",
            "per_cluster_scores": {},
            "has_substructure": {},
            "recommendation": "ok",
            "details": "Phiclust R package not installed. Run: devtools::install_github(\"semraulab/phiclust\")",
            "warning": "phiclust R package required.",
        }

    scores = {}
    has_substructure = {}

    for cluster in adata.obs[cluster_key].unique():
        mask = adata.obs[cluster_key] == cluster
        subset = adata[mask]

        if subset.n_obs < 10:
            scores[str(cluster)] = 0.0
            has_substructure[str(cluster)] = False
            continue

        try:
            # 调用 R phiclust
            expr_matrix = subset.X
            if hasattr(expr_matrix, "toarray"):
                expr_matrix = expr_matrix.toarray()

            # 构造 R 调用
            ro.globalenv["expr"] = ro.FloatMatrix(expr_matrix.T)
            ro.globalenv["clusters"] = ro.StrVector([str(cluster)] * subset.n_obs)

            result = ro.r("""
                library(phiclust)
                out <- tryCatch(
                    phiclust(expr = expr, clusters = clusters),
                    error = function(e) NULL
                )
                if (!is.null(out)) out$phi.clust else 0
            """)

            phi = float(result[0]) if result else 0.0
            scores[str(cluster)] = phi
            has_substructure[str(cluster)] = phi > threshold

        except Exception:
            scores[str(cluster)] = 0.0
            has_substructure[str(cluster)] = False

    # 综合建议
    n_sub = sum(has_substructure.values())
    total = len(has_substructure)

    if total == 0:
        recommendation = "ok"
    elif n_sub / total > 0.5:
        recommendation = "increase_resolution"
    elif n_sub == 0 and total > 5:
        recommendation = "maybe_overcluster"
    else:
        recommendation = "ok"

    return {
        "tool": "phiclust",
        "per_cluster_scores": {k: round(v, 4) for k, v in scores.items()},
        "has_substructure": has_substructure,
        "recommendation": recommendation,
        "details": f"{n_sub}/{total} clusters show hidden substructure (threshold={threshold})",
    }

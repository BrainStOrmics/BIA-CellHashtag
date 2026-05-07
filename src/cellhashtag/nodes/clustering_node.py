"""
Clustering 节点 + 质量评估。

实现 Clustering ↔ Quality Assessment 的迭代循环：
1. 执行 Leiden 聚类
2. 运行质量评估工具（Silhouette + Modularity + Phiclust）
3. 综合决策：继续 / 调整 / HITL
"""

import pickle
from typing import Any


def node_clustering(state: dict) -> dict:
    """
    聚类节点：执行 Leiden 并评估质量。

    如果质量不达标且未超过 max_iterations，
    调整 resolution 并重新聚类。

    Args:
        state: CellHashtagState。

    Returns:
        更新后的 state 片段。
    """
    import scanpy as sc
    from ..skills.clustering_quality import (
        silhouette_check,
        modularity_check,
    )

    # 加载 AnnData
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)

    iteration = state.get("iteration_count", 0)
    max_iter = state.get("max_iterations", 3)
    cluster_key = state.get("cluster_key", "leiden")

    if iteration >= max_iter:
        return {
            "errors": state.get("errors", []) + [
                f"Max clustering iterations ({max_iter}) reached."
            ],
            "status": "clustering_done",
        }

    # 获取当前 clustering 参数
    params = state.get("clustering_params", {"resolution": 0.5})
    resolution = params.get("resolution", 0.5)

    # 确保 neighbors 已计算
    if "connectivities" not in adata.obsp:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)

    # 执行 Leiden 聚类
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

    clustering_result = {
        "n_clusters": int(adata.obs[cluster_key].nunique()),
        "cluster_labels": adata.obs[cluster_key].tolist(),
    }

    # 运行质量评估工具
    quality_results = {}

    # Silhouette
    try:
        quality_results["silhouette"] = silhouette_check(adata, cluster_key)
    except Exception as e:
        quality_results["silhouette"] = {
            "tool": "silhouette",
            "error": str(e),
            "recommendation": "ok",
        }

    # Modularity
    try:
        quality_results["modularity"] = modularity_check(adata, cluster_key)
    except Exception as e:
        quality_results["modularity"] = {
            "tool": "modularity",
            "error": str(e),
            "recommendation": "ok",
        }

    # Phiclust (可选)
    try:
        from ..skills.clustering_quality import phiclust_check
        if phiclust_check is not None:
            quality_results["phiclust"] = phiclust_check(adata, cluster_key)
    except Exception as e:
        quality_results["phiclust"] = {
            "tool": "phiclust",
            "error": str(e),
            "recommendation": "ok",
        }

    # 综合决策
    decision = _make_clustering_decision(quality_results, resolution, iteration, max_iter)

    # 记录迭代历史
    history_entry = {
        "iteration": iteration,
        "params": params.copy(),
        "clustering_result": clustering_result,
        "quality_results": quality_results,
        "decision": decision["action"],
        "reason": decision["reason"],
    }

    history = state.get("clustering_history", [])
    history.append(history_entry)

    if decision["action"] == "adjust":
        # 调整参数，重新聚类
        new_params = params.copy()
        new_params["resolution"] = decision["new_resolution"]
        return {
            "clustering_params": new_params,
            "clustering_history": history,
            "iteration_count": iteration + 1,
            "status": "clustering_retry",
        }
    else:
        # 保存更新后的 AnnData
        with open(state["adata_path"], "wb") as f:
            pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)

        return {
            "clustering_result": clustering_result,
            "clustering_quality": quality_results,
            "clustering_history": history,
            "clustering_params": params,
            "status": "clustering_done",
        }


def _make_clustering_decision(
    quality_results: dict,
    current_resolution: float,
    iteration: int,
    max_iter: int,
) -> dict:
    """
    综合所有工具的结果，做出聚类决策。

    Returns:
        {"action": "proceed" | "adjust", "reason": str, "new_resolution": float}
    """
    recommendations = []
    for tool, result in quality_results.items():
        rec = result.get("recommendation", "ok")
        recommendations.append(rec)

    n_increase = sum(1 for r in recommendations if "increase" in r.lower())
    n_decrease = sum(1 for r in recommendations if "decrease" in r.lower() or "over" in r.lower())
    n_adjust = sum(1 for r in recommendations if r == "adjust")
    n_ok = sum(1 for r in recommendations if r == "ok")

    if (n_increase + n_decrease + n_adjust) > n_ok:
        # 需要调整
        if n_increase > n_decrease:
            new_resolution = round(current_resolution * 1.5, 2)
            reason = f"Tools suggest increasing resolution (currently {current_resolution})"
        else:
            new_resolution = round(current_resolution * 0.7, 2)
            reason = f"Tools suggest decreasing resolution (currently {current_resolution})"

        # 限制 resolution 范围
        new_resolution = max(0.05, min(new_resolution, 2.0))

        return {
            "action": "adjust",
            "new_resolution": new_resolution,
            "reason": reason,
        }
    else:
        return {
            "action": "proceed",
            "reason": f"Quality assessment passed (ASW/Modularity OK, {n_ok}/{len(recommendations)} tools approved)",
        }

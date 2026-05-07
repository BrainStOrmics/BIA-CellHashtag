"""
Clustering 子图。

独立的聚类 + 质量评估 LangGraph 子图，
可被主图调用，也可单独使用。

子图流程:
  START → cluster → assess_quality → [retry | finalize]
"""

from typing import TypedDict, Optional
from langgraph.graph import StateGraph, START, END


class ClusteringSubGraphState(TypedDict):
    """聚类子图状态。"""
    adata_path: str                     # AnnData 文件路径
    cluster_key: str                    # 聚类列名
    clustering_params: dict             # {resolution, n_neighbors, ...}
    clustering_result: dict             # {n_clusters, cluster_labels}
    clustering_quality: dict            # 质量评估结果
    clustering_history: list            # 迭代历史
    iteration_count: int                # 当前迭代次数
    max_iterations: int                 # 最大迭代次数
    status: str                         # clustering / assessing / done / retry
    errors: list                        # 错误列表


def build_clustering_subgraph() -> StateGraph:
    """
    构建聚类子图。

    Returns:
        编译后的 LangGraph 子图。
    """
    builder = StateGraph(ClusteringSubGraphState)

    # 节点
    builder.add_node("cluster", _node_run_clustering)
    builder.add_node("assess_quality", _node_assess_quality)

    # 边
    builder.add_edge(START, "cluster")
    builder.add_edge("cluster", "assess_quality")
    builder.add_conditional_edges(
        "assess_quality",
        _route_after_assessment,
        {
            "retry": "cluster",
            "done": END,
        },
    )

    return builder


def _node_run_clustering(state: ClusteringSubGraphState) -> dict:
    """执行 Leiden 聚类。"""
    import pickle
    import scanpy as sc

    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)

    params = state.get("clustering_params", {"resolution": 0.5})
    resolution = params.get("resolution", 0.5)
    cluster_key = state.get("cluster_key", "leiden")

    # 确保 neighbors 已计算
    if "connectivities" not in adata.obsp:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)

    # Leiden 聚类
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

    # 保存更新后的 AnnData
    with open(state["adata_path"], "wb") as f:
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)

    return {
        "clustering_result": {
            "n_clusters": int(adata.obs[cluster_key].nunique()),
            "cluster_labels": adata.obs[cluster_key].tolist(),
        },
        "iteration_count": state.get("iteration_count", 0) + 1,
    }


def _node_assess_quality(state: ClusteringSubGraphState) -> dict:
    """评估聚类质量。"""
    import pickle
    from ..skills.clustering_quality import (
        silhouette_check,
        modularity_check,
    )

    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)

    cluster_key = state.get("cluster_key", "leiden")
    quality_results = {}

    # Silhouette
    try:
        quality_results["silhouette"] = silhouette_check(adata, cluster_key)
    except Exception as e:
        quality_results["silhouette"] = {
            "tool": "silhouette", "error": str(e), "recommendation": "ok",
        }

    # Modularity
    try:
        quality_results["modularity"] = modularity_check(adata, cluster_key)
    except Exception as e:
        quality_results["modularity"] = {
            "tool": "modularity", "error": str(e), "recommendation": "ok",
        }

    # Phiclust (可选)
    try:
        from ..skills.clustering_quality import phiclust_check
        if phiclust_check is not None:
            quality_results["phiclust"] = phiclust_check(adata, cluster_key)
    except Exception:
        pass

    # 记录历史
    history = state.get("clustering_history", [])
    history.append({
        "iteration": state.get("iteration_count", 0),
        "params": state.get("clustering_params", {}),
        "quality_results": quality_results,
    })

    return {
        "clustering_quality": quality_results,
        "clustering_history": history,
    }


def _route_after_assessment(state: ClusteringSubGraphState) -> str:
    """
    根据质量评估结果路由。

    - 需要调整 → retry
    - 通过或达到最大迭代 → done
    """
    quality = state.get("clustering_quality", {})
    iteration = state.get("iteration_count", 0)
    max_iter = state.get("max_iterations", 3)

    if iteration >= max_iter:
        return "done"

    # 检查是否需要调整
    needs_adjust = False
    for tool, result in quality.items():
        rec = result.get("recommendation", "ok")
        if rec in ("adjust", "increase_resolution", "decrease_resolution", "maybe_overcluster"):
            needs_adjust = True
            break

    if needs_adjust:
        # 调整 resolution
        params = state.get("clustering_params", {}).copy()
        current_res = params.get("resolution", 0.5)

        # 根据建议调整
        phiclust = quality.get("phiclust", {})
        if phiclust.get("recommendation") == "increase_resolution":
            params["resolution"] = round(current_res * 1.5, 2)
        else:
            params["resolution"] = round(current_res * 0.7, 2)

        params["resolution"] = max(0.05, min(params["resolution"], 2.0))

        return "retry"

    return "done"

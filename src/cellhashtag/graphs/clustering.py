"""Clustering subgraph with quality assessment loop."""

import pickle
from typing import TypedDict, Optional

import scanpy as sc
from langgraph.graph import StateGraph, START, END

from ..skills.clustering_quality import silhouette_check, modularity_check


class ClusteringSubGraphState(TypedDict):
    adata_path: str
    cluster_key: str
    clustering_params: dict
    clustering_result: dict
    clustering_quality: dict
    clustering_history: list
    iteration_count: int
    max_iterations: int
    status: str
    errors: list


def build_clustering_subgraph() -> StateGraph:
    builder = StateGraph(ClusteringSubGraphState)
    builder.add_node("cluster", _run_clustering)
    builder.add_node("assess_quality", _assess_quality)
    builder.add_edge(START, "cluster")
    builder.add_edge("cluster", "assess_quality")
    builder.add_conditional_edges(
        "assess_quality", _route_after_assessment, {"retry": "cluster", "done": END}
    )
    return builder


def _run_clustering(state: ClusteringSubGraphState) -> dict:
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    cluster_key = state.get("cluster_key", "leiden")
    resolution = state.get("clustering_params", {}).get("resolution", 0.5)
    if "connectivities" not in adata.obsp:
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)
    with open(state["adata_path"], "wb") as f:
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)
    return {
        "clustering_result": {
            "n_clusters": int(adata.obs[cluster_key].nunique()),
            "cluster_labels": adata.obs[cluster_key].tolist(),
        },
        "iteration_count": state.get("iteration_count", 0) + 1,
    }


def _assess_quality(state: ClusteringSubGraphState) -> dict:
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    cluster_key = state.get("cluster_key", "leiden")
    quality = {}
    try:
        quality["silhouette"] = silhouette_check(adata, cluster_key)
    except Exception as e:
        quality["silhouette"] = {"tool": "silhouette", "error": str(e), "recommendation": "ok"}
    try:
        quality["modularity"] = modularity_check(adata, cluster_key)
    except Exception as e:
        quality["modularity"] = {"tool": "modularity", "error": str(e), "recommendation": "ok"}
    history = state.get("clustering_history", [])
    history.append({
        "iteration": state.get("iteration_count", 0),
        "params": state.get("clustering_params", {}),
        "quality": quality,
    })
    return {"clustering_quality": quality, "clustering_history": history}


def _route_after_assessment(state: ClusteringSubGraphState) -> str:
    iteration = state.get("iteration_count", 0)
    if iteration >= state.get("max_iterations", 3):
        return "done"
    quality = state.get("clustering_quality", {})
    for result in quality.values():
        rec = result.get("recommendation", "ok")
        if rec in ("adjust", "increase_resolution", "decrease_resolution", "maybe_overcluster"):
            params = state.get("clustering_params", {}).copy()
            current = params.get("resolution", 0.5)
            if quality.get("phiclust", {}).get("recommendation") == "increase_resolution":
                params["resolution"] = round(current * 1.5, 2)
            else:
                params["resolution"] = round(current * 0.7, 2)
            params["resolution"] = max(0.05, min(params["resolution"], 2.0))
            return "retry"
    return "done"

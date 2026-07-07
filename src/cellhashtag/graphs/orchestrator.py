"""Main orchestrator: perception → clustering → annotation → output."""

import os
import pickle
from pathlib import Path
from typing import Any, Optional, TypedDict, Annotated
import operator

from langgraph.graph import StateGraph, START, END

from .clustering import build_clustering_subgraph


class OrchestratorState(TypedDict):
    adata_path: str
    input_format: str
    omics_type: str
    cluster_key: str
    markers_per_cluster: Annotated[dict, operator.add]
    annotation_results: Annotated[list, operator.add]
    output_dir: str
    status: str
    errors: list


def build_orchestrator_graph() -> StateGraph:
    builder = StateGraph(OrchestratorState)
    builder.add_node("perception", node_perception)
    builder.add_node("clustering", node_clustering)
    builder.add_node("annotation", node_annotation)
    builder.add_node("output", node_output)
    builder.add_edge(START, "perception")
    builder.add_edge("perception", "clustering")
    builder.add_edge("clustering", "annotation")
    builder.add_edge("annotation", "output")
    builder.add_edge("output", END)
    return builder


def node_perception(state: OrchestratorState) -> dict[str, Any]:
    import scanpy as sc
    from ..utils.io import load_adata, perceive_data, adata_to_pickle, extract_markers

    adata = load_adata(state["adata_path"], state.get("input_format", "h5ad"))
    info = perceive_data(adata)
    if state.get("omics_type", "auto") == "auto":
        omics_type = info["omics_type"]
    else:
        omics_type = state["omics_type"]

    pickle_path = str(Path(state["adata_path"]).with_suffix(".pkl"))
    adata_to_pickle(adata, pickle_path)

    cluster_key = state.get("cluster_key", "leiden")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)

    if cluster_key not in adata.obs.columns:
        sc.tl.leiden(adata, key_added=cluster_key)

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, reference="rest")

    markers = {}
    for cluster in adata.obs[cluster_key].unique():
        markers[str(cluster)] = extract_markers(adata, cluster_key, str(cluster))

    return {
        "adata_path": pickle_path,
        "omics_type": omics_type,
        "markers_per_cluster": markers,
    }


def node_clustering(state: OrchestratorState) -> dict[str, Any]:
    import pickle
    from ..config import load_config

    config = load_config()
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    cluster_key = state.get("cluster_key", "leiden")
    subgraph = build_clustering_subgraph().compile()
    cluster_state = {
        "adata_path": state["adata_path"],
        "cluster_key": cluster_key,
        "clustering_params": {"resolution": config.clustering.default_resolution},
        "clustering_result": {},
        "clustering_quality": {},
        "clustering_history": [],
        "iteration_count": 0,
        "max_iterations": config.clustering.max_iterations,
        "status": "running",
        "errors": [],
    }
    result = subgraph.invoke(cluster_state)
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    clusters = adata.obs[cluster_key].unique().tolist()
    return {"status": f"clustering_done:{len(clusters)} clusters"}


def node_annotation(state: OrchestratorState) -> dict[str, Any]:
    import pickle
    from deepagents import create_deep_agent, SubAgent
    from ..config import load_config, setup_llm

    config = load_config()
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    cluster_key = state.get("cluster_key", "leiden")
    clusters = adata.obs[cluster_key].unique().tolist()
    markers = state.get("markers_per_cluster", {})
    skills_dir = Path(__file__).parent.parent / "skills"

    llm = setup_llm(config.llm)
    model_str = f"{config.llm.provider}:{config.llm.model}"

    deep_agent = create_deep_agent(
        model=model_str if llm is None else llm,
        system_prompt=_annotation_system_prompt(),
        subagents=[
            SubAgent(
                name="lats-search",
                description="Run MCTS tree search to find best cell type annotation for a cluster",
                system_prompt=_lats_subagent_prompt(),
                skills=[str(skills_dir / "lats_search"), str(skills_dir / "cellwiki")],
            ),
        ],
        skills=[str(skills_dir / "cellwiki"), str(skills_dir / "clustering_quality")],
    )

    results = []
    for cluster in clusters:
        cluster_id = str(cluster)
        cluster_markers = markers.get(cluster_id, [])
        try:
            agent_input = {
                "messages": [
                    {
                        "role": "user",
                        "content": _build_annotation_task(cluster_id, cluster_markers, adata),
                    }
                ]
            }
            agent_result = deep_agent.invoke(agent_input)
            results.append(_parse_agent_result(agent_result, cluster_id))
        except Exception as exc:
            results.append({
                "cluster": cluster_id,
                "cell_type": "Unknown",
                "confidence": 0.0,
                "reasoning": f"Annotation failed: {exc}",
            })

    return {
        "annotation_results": results,
        "status": "annotation_done",
    }


def node_output(state: OrchestratorState) -> dict[str, Any]:
    import pickle
    from pathlib import Path

    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)
    output_dir = Path(state.get("output_dir", "output"))
    output_dir.mkdir(parents=True, exist_ok=True)
    results = state.get("annotation_results", [])
    from ..utils.hitl import review_low_confidence_annotations, build_review_report
    results = review_low_confidence_annotations(results)
    build_review_report(results, output_dir)
    cluster_key = state.get("cluster_key", "leiden")
    for r in results:
        mask = adata.obs[cluster_key].astype(str) == r["cluster"]
        adata.obs.loc[mask, "Cell#"] = r["cell_type"]
        adata.obs.loc[mask, "Cell#_confidence"] = r["confidence"]
    out_path = output_dir / "annotated.h5ad"
    adata.write_h5ad(out_path)
    _save_csv(results, output_dir / "annotation_table.csv")
    _generate_report(results, output_dir, state.get("omics_type", "unknown"), cluster_key)
    _generate_umap(adata, output_dir, cluster_key)
    return {"status": "done", "errors": []}


_PROMPTS_DIR = Path(__file__).parent.parent / "prompts"


def load_prompt(filename: str) -> str:
    return (_PROMPTS_DIR / filename).read_text()


def _annotation_system_prompt() -> str:
    return load_prompt("annotation.md")


def _lats_subagent_prompt() -> str:
    return load_prompt("lats_search.md")


def _build_annotation_task(cluster_id: str, markers: list, adata: Any) -> str:
    import numpy as np
    expr_lines = []
    for gene in markers[:10]:
        if gene in adata.var_names:
            mat = adata[:, adata.var_names == gene].X
            if hasattr(mat, "toarray"):
                mat = mat.toarray()
            mean = float(np.mean(mat))
            expr_lines.append(f"  {gene}: mean={mean:.3f}")
    template = load_prompt("annotation_task.md")
    return template.format(
        cluster_id=cluster_id,
        markers=markers[:15],
        expression_lines="\n".join(expr_lines),
    )


def _parse_agent_result(result: dict[str, Any], cluster_id: str) -> dict[str, Any]:
    import json
    import re
    messages = result.get("messages", [])
    for msg in reversed(messages):
        content = getattr(msg, "content", "") if hasattr(msg, "content") else str(msg)
        m = re.search(r"\{.*\}", content, re.DOTALL)
        if m:
            try:
                data = json.loads(m.group())
                return {
                    "cluster": cluster_id,
                    "cell_type": data.get("cell_type", "Unknown"),
                    "confidence": float(data.get("confidence", 0.5)),
                    "reasoning": data.get("reasoning", ""),
                }
            except json.JSONDecodeError:
                pass
    return {"cluster": cluster_id, "cell_type": "Unknown", "confidence": 0.3, "reasoning": ""}


def _save_csv(results: list[dict[str, Any]], path: Path):
    import csv
    if not results:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["cluster", "cell_type", "confidence", "reasoning"])
        writer.writeheader()
        for r in results:
            writer.writerow(r)


def _generate_report(results: list[dict[str, Any]], output_dir: Path, omics_type: str, cluster_key: str):
    import datetime

    if not results:
        return

    confidences = [r["confidence"] for r in results]
    mean_conf = sum(confidences) / len(confidences)
    count_high = sum(1 for c in confidences if c >= 0.8)
    count_med = sum(1 for c in confidences if 0.5 <= c < 0.8)
    count_low = sum(1 for c in confidences if c < 0.5)

    rows = []
    for r in results:
        reasoning = r.get("reasoning", "").replace("\n", " ")
        rows.append(f"| {r['cluster']} | {r['cell_type']} | {r['confidence']:.2f} | {reasoning} |")

    md = f"""# Cell# Annotation Report

Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
Input: {omics_type} data
Clustering: {cluster_key}

## Summary

| Metric | Value |
|--------|-------|
| Total clusters annotated | {len(results)} |
| Average confidence | {mean_conf:.2f} |
| High confidence (>=0.8) | {count_high} |
| Medium confidence (0.5-0.8) | {count_med} |
| Low confidence (<0.5) | {count_low} |

## Per-Cluster Annotations

| Cluster | Cell Type | Confidence | Reasoning |
|---------|-----------|------------|-----------|
{chr(10).join(rows)}

## Recommendations

- Review low-confidence annotations manually
- Check UMAP visualization for cluster separation
"""
    (output_dir / "annotation_report.md").write_text(md)


def _generate_umap(adata: Any, output_dir: Path, cluster_key: str):
    import scanpy as sc
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            sc.pp.pca(adata)
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    sc.settings.figdir = str(output_dir)
    sc.settings.set_figure_params(dpi=150)

    sc.pl.umap(adata, color="Cell#", save="_celltype.png", show=False)
    plt.close("all")
    sc.pl.umap(adata, color="Cell#_confidence", save="_confidence.png", show=False, cmap="RdYlGn")
    plt.close("all")

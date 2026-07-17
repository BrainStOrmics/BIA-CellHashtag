"""LATS outer loop rewired onto DAG-MCTS.

Preserves run_lats_search() signature for orchestrator. Internally uses:
- cold_start (3-stage root build)
- dag_mcts (P-UCT + Soft-Bellman + confidence discount + transposition)
- tree_viz (per-step HTML dump)
- reflection + 5-component reward

DeepAgents still provides LLM calls for action_generator and eval_fn.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Optional

from .cold_start import cold_start
from .dag_mcts import DAGNode, TranspositionTable, make_state_hash, run_dag_mcts
from .tree_viz import TreeVisualizer


DEFAULT_WEIGHTS = {
    "purity": 0.25,
    "specificity": 0.20,
    "context": 0.20,
    "lats": 0.20,
    "known_marker": 0.15,
}

DEFAULT_SEARCH_PARAMS = {
    "max_iterations": 20,
    "cpuct": 1.5,
    "tau0": 1.0,
    "alpha": 0.3,
    "epsilon": 0.05,
    "early_stop_threshold": 0.95,
    "confidence_threshold": 0.7,
}


_EXPANSION_PROMPT = """\
You are generating cell type annotation actions for a single-cell RNA-seq cluster via DAG-MCTS.

Context:
- Cluster markers: {markers}
- Tissue: {tissue}
- Current node hypothesis: {current_hypothesis}
- Already-explored actions from this node: {explored_actions}

Propose 1-3 candidate actions from the DAG action space:
  - Split(res=<resolution>)  -> split current cluster further
  - Merge(target=<cluster>)  -> merge with another cluster
  - Freeze                   -> finalize current cluster label
  - AssignLabel(cell_type=<type>) -> commit to a cell type

For each action, output:
  {{"action": "<Split|Merge|Freeze|AssignLabel>", "target": "<param>", "weight": 0.0-1.0, "reasoning": "..."}}

Output ONLY valid JSON: {{"actions": [...]}}
"""

_EVALUATION_PROMPT = """\
Adversarial evaluator. Find why this hypothesis might be WRONG.

Hypothesis: {hypothesis}
Cluster markers: {markers}
Tissue: {tissue}
Depth in DAG: {depth}

Score on 5 dimensions (0.0-1.0, higher = better):
1. purity: cluster internal homogeneity (silhouette-like)
2. specificity: marker exclusivity vs other clusters
3. context: consistency with tissue ontology
4. lats: reflection confidence (self-assessment)
5. known_marker: overlap with CellWiki canonical markers

Output ONLY valid JSON:
{{"purity": 0.0, "specificity": 0.0, "context": 0.0, "lats": 0.0, "known_marker": 0.0,
  "fragmentation": 0.0,
  "reflection_output": {{"diagnosis": "...", "assessment": 0.0, "suggested_next": "...",
                        "hallucination_risk": "low|medium|high", "should_continue": true}},
  "critique_text": "..."}}
"""


def _parse_json(content: str) -> dict | None:
    m = re.search(r"\{.*\}", content, re.DOTALL)
    if m:
        try:
            return json.loads(m.group())
        except json.JSONDecodeError:
            return None
    return None


def _extract_text(result: Any) -> str:
    if isinstance(result, str):
        return result
    messages = result.get("messages", []) if isinstance(result, dict) else []
    for msg in reversed(messages):
        content = getattr(msg, "content", "") if hasattr(msg, "content") else str(msg)
        if content:
            return content
    return ""


def _tree_stats(dag: TranspositionTable) -> dict:
    nodes = dag.all_nodes()
    total = len(nodes)
    max_depth = max((n.depth for n in nodes), default=0)
    multi_parent = sum(1 for n in nodes if len(n.parents) > 1)
    return {"total_nodes": total, "max_depth": max_depth, "dag_nodes_with_multiple_parents": multi_parent}


def run_lats_search(
    cluster_id: str,
    markers: list[str],
    tissue: str = "unknown",
    skills_dir: Path | None = None,
    llm_model: str = "qwen3.5-plus",
    weights: dict[str, float] | None = None,
    search_params: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
    adata: Any | None = None,
    output_dir: str | Path = "output",
) -> dict[str, Any]:
    weights = weights or DEFAULT_WEIGHTS
    params = {**DEFAULT_SEARCH_PARAMS, **(search_params or {})}
    if skills_dir is None:
        skills_dir = Path(__file__).parent.parent / "skills"

    from deepagents import create_deep_agent

    expansion_agent = create_deep_agent(
        model=llm_model,
        system_prompt=_EXPANSION_PROMPT,
        skills=[str(skills_dir / "cellwiki"), str(skills_dir / "clustering_quality")],
    )
    eval_agent = create_deep_agent(
        model=llm_model,
        system_prompt=_EVALUATION_PROMPT,
        skills=[str(skills_dir / "cellwiki"), str(skills_dir / "clustering_quality")],
    )

    metadata = metadata or {"tissue": tissue}

    cold = cold_start(metadata, adata=adata, markers_per_cluster={cluster_id: markers})

    root_hash = make_state_hash({"cluster": cluster_id, "stage": "root"}, {"tissue": tissue})
    root = DAGNode(
        node_id=f"{cluster_id}_root",
        state_hash=root_hash,
        action=None,
        hypothesis={
            "cell_type": "UNASSIGNED",
            "expected_markers": markers[:10],
            "confidence_estimate": 0.3,
            "cluster_info": {"label": cluster_id, "heterogeneity": 0.8, "confidence": 0.3},
        },
        topology_state={"cluster": cluster_id, "stage": "root"},
        semantic_state={"tissue": tissue, "resolution_hint": cold.resolution_hint},
        depth=0,
        prior_p=1.0,
    )

    viz = TreeVisualizer(output_dir, cluster_id, enabled=params.get("enable_tree_viz", True))

    def action_generator(node: DAGNode) -> list[dict]:
        explored = []
        if node.action:
            explored.append(f"{node.action}({node.hypothesis.get('target', '')})")
        prompt_in = _EXPANSION_PROMPT.format(
            markers=", ".join(markers[:15]),
            tissue=tissue,
            current_hypothesis=json.dumps(node.hypothesis),
            explored_actions=", ".join(explored),
        )
        try:
            result = expansion_agent.invoke({"messages": [{"role": "user", "content": prompt_in}]})
            parsed = _parse_json(_extract_text(result))
            if parsed and "actions" in parsed:
                return parsed["actions"][: params.get("max_branches", 3)]
        except Exception:
            pass
        return [{"action": "AssignLabel", "target": "Unknown", "weight": 0.3}]

    def eval_fn(node: DAGNode) -> dict:
        prompt_in = _EVALUATION_PROMPT.format(
            hypothesis=json.dumps(node.hypothesis),
            markers=", ".join(markers[:15]),
            tissue=tissue,
            depth=node.depth,
        )
        try:
            result = eval_agent.invoke({"messages": [{"role": "user", "content": prompt_in}]})
            parsed = _parse_json(_extract_text(result))
            if parsed:
                return parsed
        except Exception:
            pass
        return {"purity": 0.3, "specificity": 0.3, "context": 0.3, "lats": 0.3, "known_marker": 0.3, "fragmentation": 0.1}

    def viz_callback(it: int, dag: TranspositionTable, highlight_id: str) -> None:
        viz.render_step(it, dag, highlight_id)

    result = run_dag_mcts(
        root=root,
        action_generator=action_generator,
        eval_fn=eval_fn,
        max_iterations=params["max_iterations"],
        cpuct=params["cpuct"],
        tau0=params["tau0"],
        alpha=params["alpha"],
        early_stop_threshold=params["early_stop_threshold"],
        viz_callback=viz_callback,
    )

    best_node = result["best_node"]
    best_score = result["best_score"]

    if viz.enabled and best_node is not None:
        viz.render_final(result["dag"], best_node.node_id)

    if best_score < params["confidence_threshold"] or best_node is None:
        return {
            "cluster": cluster_id,
            "cell_type": "Unknown",
            "confidence": max(best_score, 0.3),
            "reasoning": "DAG-MCTS did not reach confidence threshold",
            "evidence": [],
            "tree_stats": _tree_stats(result["dag"]),
            "cold_start": {"tissue": cold.tissue, "resolution_hint": cold.resolution_hint},
        }

    return {
        "cluster": cluster_id,
        "cell_type": best_node.hypothesis.get("target") or best_node.hypothesis.get("cell_type", "Unknown"),
        "confidence": best_score,
        "reasoning": best_node.reflection.diagnosis if best_node.reflection else "",
        "evidence": [
            {"source": "dag_mcts", "node_id": best_node.node_id, "q_value": best_node.q_value,
             "visits": best_node.visits, "depth": best_node.depth,
             "hallucination_risk": best_node.reflection.hallucination_risk if best_node.reflection else "unknown"}
        ],
        "tree_stats": _tree_stats(result["dag"]),
        "cold_start": {"tissue": cold.tissue, "resolution_hint": cold.resolution_hint},
    }

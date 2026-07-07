"""LATS outer loop: MCTS search with DeepAgents for Expansion + Evaluation.

This module wires the pure MCTS algorithm (core/mcts.py) to DeepAgents,
which provides the LLM calls for hypothesis generation (Expansion) and
adversarial scoring (Evaluation).

Architecture:
    AA node -> for _ in range(max_iterations):
                  select_node()          # pure algorithm
                  expand_with_deepagents # LLM via DeepAgents
                  evaluate_with_agents   # adversarial LLM via DeepAgents
                  backpropagate()        # pure algorithm
              -> extract_best()
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .mcts import (
    MCTSNode,
    select_node,
    backpropagate,
    extract_best,
    compute_annotation_value,
    generate_node_id,
)

# -- Default parameters --------------------------------------------------

DEFAULT_WEIGHTS = {
    "marker_overlap_hard": 0.35,
    "cellwiki_hit_score": 0.25,
    "cross_validator_score": 0.20,
    "negative_penalty": 0.10,
    "overconfidence_penalty": 0.10,
}

DEFAULT_SEARCH_PARAMS = {
    "max_iterations": 10,
    "exploration_weight": 1.414,
    "max_branches": 3,
    "early_stop_threshold": 0.95,
    "confidence_threshold": 0.70,
    "discount_gamma": 0.9,
}

# -- Prompts -------------------------------------------------------------

EXPANSION_PROMPT = """\
You are generating cell type annotation hypotheses for a single-cell RNA-seq cluster.

Context:
- Top marker genes: {markers}
- Tissue context: {tissue}
- Existing hypotheses already explored: {existing_hypotheses}

Generate 1-3 new cell type hypotheses that are:
1. Biologically plausible given the tissue and markers
2. Different from the existing hypotheses listed above
3. Each with expected marker genes and brief reasoning

Output ONLY valid JSON:
{{"hypotheses": [
  {{"cell_type": "string", "expected_markers": ["gene1", "gene2", ...], "reasoning": "string", "confidence_estimate": 0.0-1.0}}
]}}
"""

EVALUATION_PROMPT = """\
You are an adversarial evaluator. Your job is to find reasons this cell type hypothesis might be WRONG.

Hypothesis: {cell_type}
Cluster markers: {cluster_markers}
Expected markers for this type: {expected_markers}
Tissue context: {tissue}

Score on 5 dimensions (0.0-1.0 each, higher = better in ALL cases):
1. marker_overlap_hard: What fraction of expected markers appear in the cluster markers?
2. cellwiki_hit_score: Is this cell type known to exist in this tissue (per CellWiki/CellMarker)?
3. cross_validator_score: Is there diverse, independent evidence supporting this assignment?
4. negative_penalty: Absence of disconfirming markers (1 = no red flags, 0 = clear contradictions).
5. overconfidence_penalty: Calibration sanity (1 = confidence matches evidence, 0 = wildly overconfident).

Output ONLY valid JSON:
{{"dimension_scores": {{"marker_overlap_hard": 0.0, "cellwiki_hit_score": 0.0, "cross_validator_score": 0.0, "negative_penalty": 0.0, "overconfidence_penalty": 0.0}}, "critique_text": "string"}}
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
    messages = result.get("messages", []) if isinstance(result, dict) else []
    for msg in reversed(messages):
        content = getattr(msg, "content", "") if hasattr(msg, "content") else str(msg)
        if content:
            return content
    return ""


def _tree_stats(root: MCTSNode) -> dict[str, int]:
    total = 0
    max_depth = 0

    def walk(n: MCTSNode, d: int):
        nonlocal total, max_depth
        total += 1
        max_depth = max(max_depth, d)
        for c in n.children:
            walk(c, d + 1)

    walk(root, 0)
    return {"total_nodes": total, "max_depth": max_depth}


def run_lats_search(
    cluster_id: str,
    markers: list[str],
    tissue: str = "unknown",
    skills_dir: Path | None = None,
    llm_model: str = "qwen3.5-plus",
    weights: dict[str, float] | None = None,
    search_params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run LATS MCTS tree search for a single cluster.

    Returns:
        {"cell_type": str, "confidence": float, "reasoning": str,
         "evidence": list, "tree_stats": dict}
    """
    weights = weights or DEFAULT_WEIGHTS
    params = {**DEFAULT_SEARCH_PARAMS, **(search_params or {})}
    if skills_dir is None:
        skills_dir = Path(__file__).parent.parent / "skills"

    root = MCTSNode(
        node_id=generate_node_id(cluster_id, ["root"]),
        cluster_id=cluster_id,
        hypothesis={
            "cell_type": "UNASSIGNED",
            "expected_markers": markers[:10],
            "confidence_estimate": 0.3,
        },
    )

    existing_hypotheses = ["UNASSIGNED"]

    from deepagents import create_deep_agent

    expansion_prompt_tpl = EXPANSION_PROMPT
    eval_prompt_tpl = EVALUATION_PROMPT

    expansion_agent = create_deep_agent(
        model=llm_model,
        system_prompt=expansion_prompt_tpl,
        skills=[str(skills_dir / "cellwiki"), str(skills_dir / "clustering_quality")],
    )
    eval_agent = create_deep_agent(
        model=llm_model,
        system_prompt=eval_prompt_tpl,
        skills=[str(skills_dir / "cellwiki"), str(skills_dir / "clustering_quality")],
    )

    for iteration in range(params["max_iterations"]):
        node = select_node(root, params["exploration_weight"])

        if node.is_fully_expanded(params["max_branches"]):
            unvisited = [c for c in node.children if c.visits == 0]
            if unvisited:
                node = unvisited[0]
            else:
                node = max(node.children, key=lambda c: c.avg_reward)

        if not node.children or node == root:
            trajectory = node.get_trajectory_path()
            expansion_input = EXPANSION_PROMPT.format(
                markers=", ".join(markers[:15]),
                tissue=tissue,
                existing_hypotheses=", ".join(existing_hypotheses),
            )
            try:
                result = expansion_agent.invoke({
                    "messages": [{"role": "user", "content": expansion_input}]
                })
                content = _extract_text(result)
                parsed = _parse_json(content)
                if parsed and "hypotheses" in parsed:
                    for hyp in parsed["hypotheses"][: params["max_branches"]]:
                        path = [h.get("cell_type", "?") for h in trajectory] + [hyp.get("cell_type", "?")]
                        child = MCTSNode(
                            node_id=generate_node_id(cluster_id, path),
                            cluster_id=cluster_id,
                            hypothesis=hyp,
                            parent=node,
                        )
                        node.add_child(child)
                        existing_hypotheses.append(hyp.get("cell_type", "?"))
            except Exception:
                pass

        if not node.children:
            continue
        target = max(node.children, key=lambda c: c.ucb_score(params["exploration_weight"], node.visits))

        eval_input = EVALUATION_PROMPT.format(
            cell_type=target.hypothesis.get("cell_type", "Unknown"),
            cluster_markers=", ".join(markers[:15]),
            expected_markers=", ".join(target.hypothesis.get("expected_markers", [])),
            tissue=tissue,
        )
        try:
            eval_result = eval_agent.invoke({
                "messages": [{"role": "user", "content": eval_input}]
            })
            eval_content = _extract_text(eval_result)
            eval_parsed = _parse_json(eval_content)
            if eval_parsed and "dimension_scores" in eval_parsed:
                weighted_score = compute_annotation_value(eval_parsed, weights)
                target.evidence_chain.append({
                    "source": "adversarial_eval",
                    "scores": eval_parsed["dimension_scores"],
                    "weighted_score": weighted_score,
                })
                target.evaluation_result = {"weighted_score": weighted_score, **eval_parsed}
                backpropagate(target, weighted_score, params["discount_gamma"])
        except Exception:
            pass

        best_node, best_score = extract_best(root)
        if best_score >= params["early_stop_threshold"]:
            break

    best_node, best_score = extract_best(root)

    if best_score < params["confidence_threshold"] and best_node.hypothesis.get("cell_type") == "UNASSIGNED":
        return {
            "cluster": cluster_id,
            "cell_type": "Unknown",
            "confidence": max(best_score, 0.3),
            "reasoning": "LATS search did not reach confidence threshold",
            "evidence": best_node.evidence_chain,
            "tree_stats": _tree_stats(root),
        }

    return {
        "cluster": cluster_id,
        "cell_type": best_node.hypothesis.get("cell_type", "Unknown"),
        "confidence": best_score,
        "reasoning": best_node.hypothesis.get("reasoning", ""),
        "evidence": best_node.evidence_chain,
        "tree_stats": _tree_stats(root),
    }

"""DAG-MCTS Agent built on LangChain LLM.

Single compact agent object drives the full DAG-MCTS algorithm for cell type annotation.
Uses one LLM instance (via langchain init_chat_model) for propose_actions + evaluate_node.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Optional

from langchain.chat_models import init_chat_model

from ..core.cold_start import cold_start
from ..core.dag_mcts import DAGNode, TranspositionTable, make_state_hash, run_dag_mcts
from ..core.tree_viz import TreeVisualizer


_SYSTEM_PROMPT = """\
You are a DAG-MCTS agent for single-cell RNA-seq cell type annotation.
Your task: analyze cluster markers, propose annotation actions, and evaluate hypotheses.

Available actions:
- Split(res=<float>): split cluster at given resolution
- Merge(target=<cluster_id>): merge with another cluster
- Freeze: finalize current cluster label
- AssignLabel(cell_type=<type>): commit to a cell type annotation

When asked to propose actions:
- Output JSON: {"actions": [{"action": "...", "target": "...", "weight": 0.0-1.0, "reasoning": "..."}]}
- Propose 1-3 actions based on marker genes and tissue context
- Weight reflects confidence (0.0-1.0)

When asked to evaluate a hypothesis:
- Score on 5 dimensions (0.0-1.0): purity, specificity, context, lats, known_marker
- Provide reflection: diagnosis, assessment (0.0-1.0), suggested_next, hallucination_risk (low/medium/high), should_continue (bool)
- Output JSON: {"purity": 0.0, "specificity": 0.0, "context": 0.0, "lats": 0.0, "known_marker": 0.0, "fragmentation": 0.0, "reflection_output": {...}, "critique_text": "..."}

Use query_cellwiki tool to retrieve canonical marker genes for cell types.
Be adversarial: find why hypotheses might be WRONG.
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
    # Handle AIMessage (raw LLM response)
    if hasattr(result, "content"):
        content = result.content
        if isinstance(content, list):
            for block in content:
                if isinstance(block, dict) and block.get("type") == "text":
                    return block.get("text", "")
        elif content:
            return str(content)
    # Handle Anthropic multi-block response in dict format (DeepAgents)
    if isinstance(result, dict) and "messages" in result:
        for msg in reversed(result["messages"]):
            content = getattr(msg, "content", None) if hasattr(msg, "content") else msg.get("content")
            if isinstance(content, list):
                for block in content:
                    if isinstance(block, dict) and block.get("type") == "text":
                        return block.get("text", "")
            elif content:
                return str(content)
    return ""


class DAGMCTSAgent:
    """Compact DAG-MCTS agent using single DeepAgents instance."""

    def __init__(
        self,
        model: str = "qwen3.5-plus",
        model_provider: str = "openai",
        enable_tree_viz: bool = True,
        output_dir: str | Path = "output",
    ):
        self.model = model
        self.model_provider = model_provider
        self.enable_tree_viz = enable_tree_viz
        self.output_dir = Path(output_dir)

        self.llm = init_chat_model(model=model, model_provider=model_provider)

    def propose_actions(self, node: DAGNode, markers: list[str], tissue: str) -> list[dict]:
        """Invoke LLM to propose 1-3 actions for given node."""
        explored = []
        if node.action:
            explored.append(f"{node.action}({node.hypothesis.get('target', '')})")

        prompt = f"""\
{_SYSTEM_PROMPT}

Cluster markers: {', '.join(markers[:15])}
Tissue: {tissue}
Current node hypothesis: {json.dumps(node.hypothesis)}
Already-explored actions: {', '.join(explored) if explored else 'none'}

Propose 1-3 candidate actions for this node."""

        try:
            result = self.llm.invoke(prompt)
            parsed = _parse_json(_extract_text(result))
            if parsed and "actions" in parsed:
                return parsed["actions"][:3]
        except Exception:
            pass
        return [{"action": "AssignLabel", "target": "Unknown", "weight": 0.3}]

    def evaluate_node(self, node: DAGNode, markers: list[str], tissue: str) -> dict:
        """Invoke LLM to evaluate node hypothesis adversarially."""
        prompt = f"""\
{_SYSTEM_PROMPT}

Hypothesis: {json.dumps(node.hypothesis)}
Cluster markers: {', '.join(markers[:15])}
Tissue: {tissue}
Depth in DAG: {node.depth}

Evaluate this hypothesis on 5 dimensions (purity, specificity, context, lats, known_marker).
Provide reflection with hallucination_risk assessment."""

        try:
            result = self.llm.invoke(prompt)
            parsed = _parse_json(_extract_text(result))
            if parsed:
                return parsed
        except Exception:
            pass
        return {
            "purity": 0.3, "specificity": 0.3, "context": 0.3, "lats": 0.3,
            "known_marker": 0.3, "fragmentation": 0.1,
            "reflection_output": {
                "diagnosis": "evaluation failed",
                "assessment": 0.3,
                "suggested_next": "Split",
                "hallucination_risk": "high",
                "should_continue": True,
            },
        }

    def run(
        self,
        cluster_id: str,
        markers: list[str],
        tissue: str = "unknown",
        adata: Any | None = None,
        max_iterations: int = 20,
        cpuct: float = 1.5,
        tau0: float = 1.0,
        alpha: float = 0.3,
        early_stop_threshold: float = 0.95,
        confidence_threshold: float = 0.7,
        metadata: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Run full DAG-MCTS algorithm for cluster annotation."""
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

        viz = TreeVisualizer(self.output_dir, cluster_id, enabled=self.enable_tree_viz)

        def action_generator(node: DAGNode) -> list[dict]:
            return self.propose_actions(node, markers, tissue)

        def eval_fn(node: DAGNode) -> dict:
            return self.evaluate_node(node, markers, tissue)

        def viz_callback(it: int, dag: TranspositionTable, highlight_id: str) -> None:
            viz.render_step(it, dag, highlight_id)

        result = run_dag_mcts(
            root=root,
            action_generator=action_generator,
            eval_fn=eval_fn,
            max_iterations=max_iterations,
            cpuct=cpuct,
            tau0=tau0,
            alpha=alpha,
            early_stop_threshold=early_stop_threshold,
            viz_callback=viz_callback,
        )

        best_node = result["best_node"]
        best_score = result["best_score"]

        if viz.enabled and best_node is not None:
            viz.render_final(result["dag"], best_node.node_id)

        if best_score < confidence_threshold or best_node is None:
            return {
                "cluster": cluster_id,
                "cell_type": "Unknown",
                "confidence": max(best_score, 0.3),
                "reasoning": "DAG-MCTS did not reach confidence threshold",
                "evidence": [],
                "tree_stats": self._tree_stats(result["dag"]),
                "cold_start": {"tissue": cold.tissue, "resolution_hint": cold.resolution_hint},
            }

        return {
            "cluster": cluster_id,
            "cell_type": best_node.hypothesis.get("target") or best_node.hypothesis.get("cell_type", "Unknown"),
            "confidence": best_score,
            "reasoning": best_node.reflection.diagnosis if best_node.reflection else "",
            "evidence": [
                {
                    "source": "dag_mcts",
                    "node_id": best_node.node_id,
                    "q_value": best_node.q_value,
                    "visits": best_node.visits,
                    "depth": best_node.depth,
                    "hallucination_risk": best_node.reflection.hallucination_risk if best_node.reflection else "unknown",
                }
            ],
            "tree_stats": self._tree_stats(result["dag"]),
            "cold_start": {"tissue": cold.tissue, "resolution_hint": cold.resolution_hint},
        }

    @staticmethod
    def _tree_stats(dag: TranspositionTable) -> dict:
        nodes = dag.all_nodes()
        total = len(nodes)
        max_depth = max((n.depth for n in nodes), default=0)
        multi_parent = sum(1 for n in nodes if len(n.parents) > 1)
        return {"total_nodes": total, "max_depth": max_depth, "dag_nodes_with_multiple_parents": multi_parent}

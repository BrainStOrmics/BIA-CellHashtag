"""Tree visualizer: DAG -> NetworkX -> PyVis HTML per MCTS iteration.

Emits interactive HTML file per iteration. Node color = hallucination_risk,
size = log(visits). Edge width = prior_P.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Optional


RISK_COLORS = {"low": "#4CAF50", "medium": "#FFC107", "high": "#F44336"}
UNVISITED_COLOR = "#9E9E9E"
HIGHLIGHT_COLOR = "#2196F3"


class TreeVisualizer:
    def __init__(self, output_dir: str | Path, cluster_id: str, enabled: bool = True):
        self.output_dir = Path(output_dir) / "trees" / cluster_id
        self.cluster_id = cluster_id
        self.enabled = enabled
        if self.enabled:
            self.output_dir.mkdir(parents=True, exist_ok=True)

    def render_step(
        self,
        iteration: int,
        dag: Any,
        highlight_id: Optional[str] = None,
    ) -> Optional[Path]:
        if not self.enabled:
            return None
        try:
            import networkx as nx
            from pyvis.network import Network
        except ImportError as e:
            print(f"tree_viz: missing dependency {e}. Install with: pip install networkx pyvis")
            self.enabled = False
            return None

        G = nx.DiGraph()
        for node in dag.all_nodes():
            risk = node.reflection.hallucination_risk if node.reflection else "low"
            color = RISK_COLORS.get(risk, UNVISITED_COLOR if node.visits == 0 else "#607D8B")
            if highlight_id and node.node_id == highlight_id:
                color = HIGHLIGHT_COLOR
            size = 10 + 8 * math.log1p(node.visits)
            label = _short_label(node)
            title = (
                f"node_id: {node.node_id}\n"
                f"action: {node.action}\n"
                f"visits: {node.visits}\n"
                f"q_value: {node.q_value:.3f}\n"
                f"prior_p: {node.prior_p:.3f}\n"
                f"hallucination_risk: {risk}\n"
                f"depth: {node.depth}\n"
                f"frozen: {node.frozen}\n"
                f"diagnosis: {(node.reflection.diagnosis if node.reflection else '')[:200]}"
            )
            G.add_node(
                node.node_id,
                label=label,
                title=title,
                color=color,
                size=size,
                shape="dot",
                font={"size": 10},
            )
            for child_hash in node.children:
                child = dag.lookup(child_hash)
                if child is not None:
                    width = 0.5 + 3.0 * child.prior_p
                    G.add_edge(
                        node.node_id,
                        child.node_id,
                        title=f"prior_P={child.prior_p:.3f}",
                        width=width,
                        color={"color": "#848484"},
                    )

        net = Network(height="700px", width="100%", directed=True, notebook=False)
        net.from_nx(G)
        net.set_options(_VIS_OPTIONS)

        out_path = self.output_dir / f"iter_{iteration:03d}.html"
        net.write_html(str(out_path), notebook=False, open_browser=False)
        return out_path

    def render_final(self, dag: Any, best_id: Optional[str]) -> Optional[Path]:
        return self.render_step(9999, dag, highlight_id=best_id)


def _short_label(node: Any) -> str:
    hyp = node.hypothesis or {}
    action = node.action or "?"
    target = hyp.get("target") or hyp.get("cell_type") or hyp.get("label") or ""
    if target:
        return f"{action}\n{target}"[:30]
    return f"{action}"[:20]


_VIS_OPTIONS = """
{
  "nodes": {
    "font": {"face": "monospace", "align": "center"},
    "borderWidth": 2,
    "borderWidthSelected": 3
  },
  "edges": {
    "smooth": {"type": "dynamic"},
    "arrows": {"to": {"enabled": true, "scaleFactor": 0.5}}
  },
  "physics": {
    "enabled": true,
    "hierarchicalRepulsion": {
      "centralGravity": 0.0,
      "springLength": 120,
      "springConstant": 0.01,
      "nodeDistance": 140,
      "damping": 0.09
    },
    "solver": "hierarchicalRepulsion"
  },
  "layout": {"hierarchical": {"enabled": true, "direction": "UD", "sortMethod": "directed"}}
}
"""

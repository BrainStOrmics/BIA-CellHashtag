"""Pure MCTS algorithm for LATS-based cell annotation. No LLM, no I/O."""

import math
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class MCTSNode:
    node_id: str
    cluster_id: str
    hypothesis: dict
    parent: Optional["MCTSNode"] = None
    children: list["MCTSNode"] = field(default_factory=list)
    visits: int = 0
    total_reward: float = 0.0
    evidence_chain: list[dict] = field(default_factory=list)
    evaluation_result: Optional[dict] = None

    @property
    def avg_reward(self) -> float:
        return self.total_reward / self.visits if self.visits > 0 else 0.0

    @property
    def confidence(self) -> float:
        if self.evaluation_result and "weighted_score" in self.evaluation_result:
            return self.evaluation_result["weighted_score"]
        return self.hypothesis.get("confidence_estimate", 0.5)

    def ucb_score(self, exploration_weight: float = 1.414, parent_visits: int | None = None) -> float:
        if parent_visits is None:
            parent_visits = self.parent.visits if self.parent else 1

        if self.visits == 0:
            return float("inf")

        exploitation = self.avg_reward
        exploration = exploration_weight * math.sqrt(math.log(parent_visits) / self.visits)

        return exploitation + exploration

    def add_child(self, child: "MCTSNode"):
        child.parent = self
        self.children.append(child)

    def is_fully_expanded(self, max_branches: int) -> bool:
        return len(self.children) >= max_branches

    def get_trajectory_path(self) -> list[dict]:
        """Return hypothesis path from root to this node for DeepAgents context."""
        path = []
        node = self
        while node:
            path.append(node.hypothesis)
            node = node.parent
        return path[::-1]


def select_node(root: MCTSNode, exploration_weight: float = 1.414, max_branches: int = 3) -> MCTSNode:
    """Select a leaf node for expansion.

    Walks down choosing the highest-UCT child, but stops at the first node
    that is NOT fully expanded (has room for more children). This prevents
    the tree from growing into long chains when a flat hypothesis space is intended.
    """
    current = root
    while current.children and current.is_fully_expanded(max_branches):
        current = max(current.children, key=lambda c: c.ucb_score(exploration_weight, current.visits))
    unvisited = [c for c in current.children if c.visits == 0]
    if unvisited:
        return max(unvisited, key=lambda c: c.ucb_score(exploration_weight, current.visits))
    return current


def backpropagate(node: MCTSNode, reward: float, gamma: float = 1.0):
    """Propagate reward up the tree with optional discount for depth."""
    current = node
    depth = 0
    while current:
        current.visits += 1
        current.total_reward += reward * (gamma ** depth)
        current = current.parent
        depth += 1


def extract_best(root: MCTSNode) -> tuple[MCTSNode, float]:
    best = root

    def traverse(n: MCTSNode):
        nonlocal best
        if n.avg_reward > best.avg_reward:
            best = n
        for c in n.children:
            traverse(c)

    traverse(root)
    return best, best.avg_reward


def compute_annotation_value(evaluation: dict, weights: dict[str, float]) -> float:
    dims = evaluation.get("dimension_scores", {})
    return sum(dims.get(dim, 0.0) * w for dim, w in weights.items())


def generate_node_id(cluster_id: str, path: list[str]) -> str:
    path_hash = hash(tuple(path)) & 0xFFFF
    return f"{cluster_id}_n{path_hash:04x}"

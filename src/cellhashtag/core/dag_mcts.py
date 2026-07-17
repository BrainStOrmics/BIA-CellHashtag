"""DAG-MCTS core: P-UCT selection, Soft-Bellman backprop, confidence discounting, transposition.

Replaces the tree-only core/mcts.py. Search space is a DAG — cross-branch Merge
and ReEmbed create multi-parent nodes. Transposition table dedups convergent states.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from .reflection import ReflectionSchema, parse_reflection
from .reward import (
    ClusterStats,
    RewardBreakdown,
    SlidingWindowNormalizer,
    compute_reward,
    final_reward,
    potential,
)


@dataclass
class DAGNode:
    node_id: str
    state_hash: str
    action: Optional[str] = None
    hypothesis: dict = field(default_factory=dict)
    parents: list[str] = field(default_factory=list)
    children: list[str] = field(default_factory=list)
    visits: int = 0
    q_value: float = 0.0
    prior_p: float = 0.0
    in_flight: int = 0
    reflection: Optional[ReflectionSchema] = None
    reward_breakdown: Optional[RewardBreakdown] = None
    topology_state: dict = field(default_factory=dict)
    semantic_state: dict = field(default_factory=dict)
    depth: int = 0
    frozen: bool = False

    @property
    def avg_reward(self) -> float:
        return self.q_value


class TranspositionTable:
    def __init__(self):
        self._table: dict[str, DAGNode] = {}

    def lookup(self, state_hash: str) -> Optional[DAGNode]:
        return self._table.get(state_hash)

    def store(self, node: DAGNode) -> None:
        self._table[node.state_hash] = node

    def add_parent(self, state_hash: str, parent_id: str) -> None:
        node = self._table.get(state_hash)
        if node is not None and parent_id not in node.parents:
            node.parents.append(parent_id)

    def __len__(self) -> int:
        return len(self._table)

    def all_nodes(self) -> list[DAGNode]:
        return list(self._table.values())


def make_state_hash(topology: dict, semantic: dict) -> str:
    blob = f"{sorted(topology.items())}|{sorted(semantic.items())}".encode()
    return hashlib.sha256(blob).hexdigest()[:32]


def priority(cluster_info: dict, research_context: Optional[dict] = None) -> float:
    relevance = 1.0
    if research_context:
        targets = research_context.get("target_lineages", [])
        label = cluster_info.get("label", "").lower()
        if targets and any(t.lower() in label for t in targets):
            relevance = 1.0
        elif targets:
            relevance = 0.3
    heterogeneity = float(cluster_info.get("heterogeneity", 0.5))
    uncertainty = 1.0 - float(cluster_info.get("confidence", 0.5))
    return relevance * heterogeneity * uncertainty


def fpu(node: DAGNode, v_root: float, delta: float, priority_score: float, theta_high: float, theta_low: float) -> float:
    if node.visits > 0:
        return node.q_value
    if priority_score >= theta_high:
        return v_root + delta
    if priority_score < theta_low:
        return v_root - delta
    return v_root


def select_node(
    dag: TranspositionTable,
    root_id: str,
    cpuct: float = 1.5,
    research_context: Optional[dict] = None,
    v_root: float = 0.5,
    fpu_delta: float = 0.5,
    theta_high: float = 0.7,
    theta_low: float = 0.3,
) -> DAGNode:
    current = dag.lookup(root_id)
    if current is None:
        raise ValueError(f"Root {root_id} missing")
    seen: set[str] = {current.node_id}
    while current.children and all(dag.lookup(c) is not None for c in current.children):
        child_nodes = [dag.lookup(c) for c in current.children if dag.lookup(c) is not None]
        if not child_nodes:
            break
        scored = []
        for c in child_nodes:
            p = priority(c.hypothesis.get("cluster_info", {}), research_context)
            urg = fpu(c, v_root, fpu_delta, p, theta_high, theta_low)
            ucb = urg + cpuct * c.prior_p * math.sqrt(max(current.visits, 1)) / (1 + c.visits)
            scored.append((ucb, c))
        current = max(scored, key=lambda x: x[0])[1]
        if current.frozen or current.node_id in seen:
            break
        seen.add(current.node_id)
    return current


def expand(
    node: DAGNode,
    dag: TranspositionTable,
    action_shortlist: list[dict],
    epsilon: float = 0.05,
) -> list[DAGNode]:
    new_children: list[DAGNode] = []
    total_conf = sum(a.get("weight", 0.5) for a in action_shortlist) or 1.0
    for act in action_shortlist:
        prior = (act.get("weight", 0.5) / total_conf) * (1 - epsilon) + epsilon / max(len(action_shortlist), 1)
        child_state_hash = make_state_hash(
            {**node.topology_state, "action": act.get("action", "?")},
            {**node.semantic_state, "target": act.get("target", "?")},
        )
        if child_state_hash == node.state_hash:
            continue
        existing = dag.lookup(child_state_hash)
        if existing is not None:
            dag.add_parent(child_state_hash, node.state_hash)
            if node.state_hash not in existing.parents:
                existing.parents.append(node.state_hash)
            if child_state_hash not in node.children:
                node.children.append(child_state_hash)
            new_children.append(existing)
            continue
        child = DAGNode(
            node_id=f"{node.node_id}->{act.get('action', '?')}:{act.get('target', '?')}",
            state_hash=child_state_hash,
            action=act.get("action"),
            hypothesis={**node.hypothesis, "action": act.get("action"), "target": act.get("target")},
            parents=[node.state_hash],
            prior_p=prior,
            depth=node.depth + 1,
            topology_state={**node.topology_state, "action": act.get("action")},
            semantic_state={**node.semantic_state, "target": act.get("target")},
        )
        dag.store(child)
        if child_state_hash not in node.children:
            node.children.append(child_state_hash)
        new_children.append(child)
    return new_children


def evaluate(
    node: DAGNode,
    eval_fn: Callable[[DAGNode], dict],
    normalizer: SlidingWindowNormalizer,
    prev_phi: float,
    weights: Optional[dict] = None,
) -> RewardBreakdown:
    scores = eval_fn(node)
    breakdown = compute_reward(
        purity=float(scores.get("purity", 0.0)),
        specificity=float(scores.get("specificity", 0.0)),
        context_score=float(scores.get("context", 0.0)),
        lats_score=float(scores.get("lats", 0.0)),
        known_marker_score=float(scores.get("known_marker", 0.0)),
        overfragmentation_penalty=float(scores.get("fragmentation", 0.0)),
        weights=weights,
    )
    normalizer.push(breakdown.raw)
    normalized = normalizer.normalize(breakdown.raw)
    new_clusters = [ClusterStats(label=str(i), silhouette=0.5, marker_coherence=0.5) for i in range(max(1, len(node.children) or 1))]
    next_phi = potential(new_clusters)
    final_reward(breakdown, normalized, next_phi, prev_phi)

    reflection = parse_reflection(scores.get("reflection_output"))
    if reflection is not None:
        node.reflection = reflection

    node.reward_breakdown = breakdown
    return breakdown


def temperature_at_depth(depth: int, tau0: float = 1.0, gamma: float = 0.25) -> float:
    return tau0 / (1.0 + gamma * depth)


def soft_bellman_v(node: DAGNode, dag: TranspositionTable, tau: float = 1.0) -> float:
    if not node.children:
        return node.q_value
    child_nodes = [dag.lookup(c) for c in node.children if dag.lookup(c) is not None]
    if not child_nodes:
        return node.q_value
    n_eff = sum(max(cn.visits, 1) for cn in child_nodes)
    log_sum = None
    for cn in child_nodes:
        w = (max(cn.visits, 1) + 1e-8) / (n_eff + 1e-8 * len(child_nodes))
        term = w * math.exp(cn.q_value / max(tau, 1e-6))
        log_sum = term if log_sum is None else log_sum + term
    if log_sum is None or log_sum <= 0:
        return node.q_value
    return tau * math.log(log_sum)


def shared_discount(node: DAGNode, alpha: float = 0.3) -> float:
    k = max(len(node.parents), 1)
    return 1.0 / (1.0 + alpha * (k - 1))


def backpropagate(
    leaf: DAGNode,
    dag: TranspositionTable,
    immediate_reward: float,
    tau0: float = 1.0,
    alpha: float = 0.3,
) -> None:
    path: list[DAGNode] = []
    current: Optional[DAGNode] = leaf
    visited_ids: set[str] = set()
    while current is not None and current.node_id not in visited_ids:
        path.append(current)
        visited_ids.add(current.node_id)
        if not current.parents:
            break
        current = dag.lookup(current.parents[0]) if current.parents else None

    for depth, node in enumerate(path):
        tau = temperature_at_depth(depth, tau0)
        v_successor = soft_bellman_v(node, dag, tau) if depth > 0 else immediate_reward
        lam = shared_discount(node, alpha) * (node.reflection.risk_weight() if node.reflection else 1.0)
        delta = immediate_reward + lam * v_successor - node.q_value
        node.visits += 1
        node.q_value = node.q_value + delta / max(node.visits, 1)


def extract_best(dag: TranspositionTable, root_id: str) -> tuple[Optional[DAGNode], float]:
    best_node: Optional[DAGNode] = None
    best_score = -float("inf")
    for node in dag.all_nodes():
        if node.visits == 0:
            continue
        if node.q_value > best_score:
            best_score = node.q_value
            best_node = node
    return best_node, best_score


def run_dag_mcts(
    *,
    root: DAGNode,
    action_generator: Callable[[DAGNode], list[dict]],
    eval_fn: Callable[[DAGNode], dict],
    max_iterations: int = 20,
    cpuct: float = 1.5,
    tau0: float = 1.0,
    alpha: float = 0.3,
    early_stop_threshold: float = 0.95,
    research_context: Optional[dict] = None,
    viz_callback: Optional[Callable[[int, TranspositionTable, str], None]] = None,
) -> dict:
    dag = TranspositionTable()
    dag.store(root)
    normalizer = SlidingWindowNormalizer()
    prev_phi = potential([ClusterStats(label="root", silhouette=0.3, marker_coherence=0.3)])

    stats = {"iterations": 0, "dag_size": 1}

    for it in range(max_iterations):
        selected = select_node(dag, root.state_hash, cpuct, research_context)

        actions = action_generator(selected)
        if actions:
            expand(selected, dag, actions)

        targets = [dag.lookup(c) for c in selected.children
                   if dag.lookup(c) is not None and dag.lookup(c).visits == 0]
        if not targets:
            if selected.visits > 0:
                continue
            targets = [selected]

        for target in targets[:1]:
            target.in_flight += 1
            breakdown = evaluate(target, eval_fn, normalizer, prev_phi)
            prev_phi = potential([ClusterStats(label="cur", silhouette=breakdown.purity, marker_coherence=breakdown.known_marker)])
            backpropagate(target, dag, breakdown.final, tau0, alpha)
            target.in_flight -= 1

        if viz_callback is not None:
            try:
                viz_callback(it, dag, selected.node_id)
            except Exception:
                pass

        stats["iterations"] = it + 1
        stats["dag_size"] = len(dag)

        best_node, best_score = extract_best(dag, root.node_id)
        if best_score >= early_stop_threshold:
            break

    best_node, best_score = extract_best(dag, root.node_id)
    return {
        "best_node": best_node,
        "best_score": best_score,
        "dag": dag,
        "stats": stats,
        "normalizer_stats": normalizer.stats(),
    }

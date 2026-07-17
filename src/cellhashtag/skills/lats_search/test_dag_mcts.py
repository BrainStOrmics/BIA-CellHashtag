"""Unit tests for DAG-MCTS core. Run with: pytest .../test_dag_mcts.py -v"""

import pytest

from cellhashtag.core.dag_mcts import (
    DAGNode,
    TranspositionTable,
    backpropagate,
    expand,
    extract_best,
    fpu,
    make_state_hash,
    priority,
    run_dag_mcts,
    select_node,
    shared_discount,
    soft_bellman_v,
    temperature_at_depth,
)
from cellhashtag.core.reflection import ReflectionSchema


def _node(node_id="n1", **kwargs):
    defaults = dict(node_id=node_id, state_hash=f"hash_{node_id}", action=None, hypothesis={}, prior_p=0.5)
    defaults.update(kwargs)
    return DAGNode(**defaults)


def test_dagnode_defaults():
    n = _node()
    assert n.visits == 0
    assert n.q_value == 0.0
    assert n.parents == []
    assert n.children == []


def test_transposition_table_store_lookup():
    tt = TranspositionTable()
    n = _node("a", state_hash="h_a")
    tt.store(n)
    assert tt.lookup("h_a") is n
    assert tt.lookup("h_missing") is None
    assert len(tt) == 1


def test_transposition_add_parent_dedups():
    tt = TranspositionTable()
    n = _node("a", state_hash="h_a")
    tt.store(n)
    tt.add_parent("h_a", "p1")
    tt.add_parent("h_a", "p1")
    tt.add_parent("h_a", "p2")
    assert n.parents == ["p1", "p2"]


def test_make_state_hash_deterministic_and_differs():
    h1 = make_state_hash({"cluster": "c1"}, {"tissue": "pbmc"})
    h2 = make_state_hash({"cluster": "c1"}, {"tissue": "pbmc"})
    h3 = make_state_hash({"cluster": "c2"}, {"tissue": "pbmc"})
    assert h1 == h2
    assert h1 != h3


def test_priority_formula():
    cluster = {"label": "T cell progenitor", "heterogeneity": 0.8, "confidence": 0.4}
    ctx = {"target_lineages": ["t cell"]}
    p = priority(cluster, ctx)
    assert abs(p - (1.0 * 0.8 * 0.6)) < 1e-6


def test_fpu_unvisited_priority_levels():
    node = _node("u", visits=0, q_value=0.0)
    assert fpu(node, 0.5, delta=0.5, priority_score=0.8, theta_high=0.7, theta_low=0.3) == 1.0
    assert fpu(node, 0.5, delta=0.5, priority_score=0.2, theta_high=0.7, theta_low=0.3) == 0.0
    assert fpu(node, 0.5, delta=0.5, priority_score=0.5, theta_high=0.7, theta_low=0.3) == 0.5


def test_fpu_visited_returns_q_value():
    node = _node("v", visits=3, q_value=0.7)
    assert fpu(node, 0.5, delta=0.5, priority_score=0.9, theta_high=0.7, theta_low=0.3) == 0.7


def test_temperature_at_depth_decays():
    t0 = temperature_at_depth(0, tau0=1.0, gamma=0.25)
    t5 = temperature_at_depth(5, tau0=1.0, gamma=0.25)
    assert t0 == 1.0
    assert abs(t5 - 1.0 / (1.0 + 0.25 * 5)) < 1e-9


def test_shared_discount_single_parent():
    n = _node("x", parents=["p1"])
    assert shared_discount(n, alpha=0.3) == 1.0


def test_shared_discount_multi_parent():
    n = _node("x", parents=["p1", "p2", "p3"])
    assert abs(shared_discount(n, alpha=0.3) - (1.0 / (1.0 + 0.3 * 2))) < 1e-9


def test_expand_creates_children_and_registers_transposition():
    tt = TranspositionTable()
    root = _node("root", state_hash="h_root")
    tt.store(root)
    actions = [{"action": "Split", "target": "0.5", "weight": 0.6}, {"action": "AssignLabel", "target": "Tcell", "weight": 0.4}]
    children = expand(root, tt, actions, epsilon=0.05)
    assert len(children) == 2
    assert len(root.children) == 2
    assert all(c.parents == ["h_root"] for c in children)
    priors = [c.prior_p for c in children]
    assert sum(priors) <= 1.0 + 1e-9


def test_expand_transposition_reuses_existing_state():
    tt = TranspositionTable()
    root1 = _node("r1", state_hash="h_r1")
    root2 = _node("r2", state_hash="h_r2")
    tt.store(root1); tt.store(root2)
    action = [{"action": "Split", "target": "0.5", "weight": 0.5}]
    c1 = expand(root1, tt, action)[0]
    c2 = expand(root2, tt, action)[0]
    assert c1 is c2
    assert len(c1.parents) == 2
    assert "h_r2" in c1.parents


def test_select_node_puct_walks_to_leaf():
    tt = TranspositionTable()
    root = _node("root", state_hash="h_root", visits=10)
    child = _node("c1", state_hash="h_c1", visits=5, q_value=0.6, prior_p=0.8)
    root.children = ["h_c1"]
    child.parents = ["h_root"]
    tt.store(root); tt.store(child)
    selected = select_node(tt, "h_root", cpuct=1.5)
    assert selected is child


def test_soft_bellman_leaf_returns_q():
    tt = TranspositionTable()
    leaf = _node("L", state_hash="h_L", q_value=0.7)
    tt.store(leaf)
    assert soft_bellman_v(leaf, tt, tau=1.0) == 0.7


def test_backpropagate_increments_visits():
    tt = TranspositionTable()
    root = _node("root", state_hash="h_root")
    child = _node("c1", state_hash="h_c1", parents=["h_root"])
    root.children = ["h_c1"]
    tt.store(root); tt.store(child)
    backpropagate(child, tt, immediate_reward=0.8, tau0=1.0, alpha=0.3)
    assert child.visits == 1
    assert root.visits == 1


def test_extract_best_picks_highest_q():
    tt = TranspositionTable()
    root = _node("root", state_hash="h_root", visits=1, q_value=0.1)
    tt.store(root)
    for i, q in enumerate([0.3, 0.9, 0.5]):
        n = _node(f"c{i}", state_hash=f"h_c{i}", visits=1, q_value=q)
        tt.store(n)
    best, score = extract_best(tt, "root")
    assert best.node_id == "c1"
    assert abs(score - 0.9) < 1e-9


def test_run_dag_mcts_end_to_end_with_dummy_callbacks():
    root = _node("root", state_hash="h_root", prior_p=1.0)
    viz_calls = []

    def action_gen(node):
        return [{"action": "Split", "target": f"d{node.depth}", "weight": 0.5}]

    def eval_fn(node):
        return {"purity": 0.6, "specificity": 0.5, "context": 0.5, "lats": 0.5, "known_marker": 0.5, "fragmentation": 0.1}

    result = run_dag_mcts(
        root=root,
        action_generator=action_gen,
        eval_fn=eval_fn,
        max_iterations=3,
        cpuct=1.5,
        tau0=1.0,
        alpha=0.3,
        early_stop_threshold=2.0,
        viz_callback=lambda it, dag, h: viz_calls.append(it),
    )
    assert result["stats"]["iterations"] == 3
    assert result["dag"].all_nodes()
    assert len(viz_calls) == 3
    assert result["best_node"] is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

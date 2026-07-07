"""MCTS scorers — pure algorithmic comparisons."""
from __future__ import annotations

import math


def best_action_match(fixture, actual: dict) -> float:
    if actual.get("best_action") == fixture.expected["best_action"]:
        return 1.0
    if actual.get("best_action") == fixture.expected.get("second_best_action"):
        return 0.40
    return 0.0


def depth_effectiveness(fixture, actual: dict) -> float:
    depth = actual.get("best_node_depth", 0)
    max_eff = fixture.expected["max_depth_effective"]
    if depth <= max_eff:
        return 1.0
    if depth == max_eff + 1:
        return 0.70
    return 0.30


def exploration_share(fixture, actual: dict) -> float:
    share = float(actual.get("exploration_share", 0.0))
    expected_min = float(fixture.expected["exploration_share_min"])
    if share >= expected_min:
        return 1.0
    if share >= expected_min - 0.10:
        return 0.70
    if share < expected_min - 0.20:
        return 0.20
    return 0.40


def ucb_consistency(fixture, actual: dict) -> float:
    """Verifies select_node matches closed-form UCB1."""
    nodes = {n["node_id"]: n for n in fixture.input["nodes"]}
    root = nodes[fixture.input["root"]["node_id"]]
    c = fixture.input["ucb_c"]
    N = root["visits"] or 1

    def ucb(child_id: str) -> float:
        child = nodes[child_id]
        if child["visits"] == 0:
            return math.inf
        exploit = child["reward"] / child["visits"]
        explore = c * math.sqrt(math.log(N) / child["visits"])
        return exploit + explore

    expected = max(root["children"], key=ucb)
    return 1.0 if actual.get("selected") == expected else 0.0


def backprop_discount(fixture, actual: dict) -> float:
    expected_root_reward = actual.get("expected_root_reward")
    computed = actual.get("computed_root_reward")
    if expected_root_reward is None or computed is None:
        return 0.0
    return 1.0 if abs(expected_root_reward - computed) < 1e-6 else 0.0


def top2_gap(fixture, actual: dict) -> float:
    values = sorted(actual.get("action_values", {}).values(), reverse=True)
    if len(values) < 2:
        return 1.0
    return values[0] - values[1]

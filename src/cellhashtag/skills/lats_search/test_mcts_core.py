"""Unit tests for MCTS core utilities. Run with: pytest .../test_mcts_core.py -v"""

import pytest

from cellhashtag.core.mcts import (
    MCTSNode,
    select_node,
    backpropagate,
    extract_best,
    compute_annotation_value,
    generate_node_id,
)


HYP = {"cell_type": "T cell", "expected_markers": ["CD3D", "CD4"], "confidence_estimate": 0.75}


class TestMCTSNode:
    def test_initial_stats(self):
        node = MCTSNode("n1", "cluster_001", HYP)
        assert node.visits == 0
        assert node.avg_reward == 0.0

    def test_update_stats(self):
        node = MCTSNode("n1", "cluster_001", HYP)
        backpropagate(node, 0.8)
        node.evaluation_result = {"weighted_score": 0.8}
        assert node.visits == 1
        assert node.avg_reward == 0.8
        assert node.confidence == 0.8

    def test_ucb_unvisited(self):
        parent = MCTSNode("root", "c1", HYP, visits=10)
        child = MCTSNode("c1", "c1", HYP, parent=parent)
        assert child.ucb_score() == float("inf")

    def test_ucb_exploitation_exploration(self):
        parent = MCTSNode("root", "c1", HYP, visits=100)
        child1 = MCTSNode("c1", "c1", HYP, parent=parent, visits=50, total_reward=40.0)
        child2 = MCTSNode("c2", "c1", HYP, parent=parent, visits=30, total_reward=10.0)
        assert child1.ucb_score(1.414, 100) > child2.ucb_score(1.414, 100)

    def test_backpropagate(self):
        root = MCTSNode("root", "c1", HYP)
        child = MCTSNode("c1", "c1", HYP)
        root.add_child(child)
        backpropagate(child, 0.9)
        assert child.visits == 1
        assert child.total_reward == 0.9
        assert root.visits == 1
        assert root.total_reward == 0.9


class TestSelectionAndExtraction:
    def test_select_node_prefers_unvisited(self):
        root = MCTSNode("root", "c1", HYP, visits=10)
        visited = MCTSNode("v1", "c1", HYP, parent=root, visits=5, total_reward=4.0)
        unvisited = MCTSNode("u1", "c1", HYP, parent=root)
        root.children = [visited, unvisited]
        assert select_node(root) == unvisited

    def test_extract_best(self):
        root = MCTSNode("root", "c1", HYP)
        for reward in [0.3, 0.9, 0.5]:
            child = MCTSNode(f"c{reward}", "c1", HYP, parent=root)
            backpropagate(child, reward)
            child.evaluation_result = {"weighted_score": reward}
            root.add_child(child)
        best, score = extract_best(root)
        assert best.avg_reward == 0.9


class TestUtilityFunctions:
    def test_compute_annotation_value(self):
        evaluation = {
            "dimension_scores": {
                "marker_match": 0.9,
                "ontology_consistency": 0.8,
                "evidence_diversity": 0.7,
                "llm_confidence": 0.6,
            }
        }
        weights = {"marker_match": 0.4, "ontology_consistency": 0.3, "evidence_diversity": 0.2, "llm_confidence": 0.1}
        score = compute_annotation_value(evaluation, weights)
        assert abs(score - (0.4 * 0.9 + 0.3 * 0.8 + 0.2 * 0.7 + 0.1 * 0.6)) < 0.001

    def test_generate_node_id_deterministic(self):
        id1 = generate_node_id("cluster_001", ["root", "T_cell"])
        id2 = generate_node_id("cluster_001", ["root", "T_cell"])
        assert id1 == id2

    def test_generate_node_id_unique(self):
        id1 = generate_node_id("cluster_001", ["root", "T_cell"])
        id2 = generate_node_id("cluster_001", ["root", "B_cell"])
        assert id1 != id2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

"""
Unit tests for MCTS core utilities.

Run with: python -m pytest skills/lats_search/test_mcts_core.py -v
"""
import pytest
from .mcts_core import (
    MCTSNode, AnnotationHypothesis, EvidenceItem,
    select_node, backpropagate, extract_best,
    compute_annotation_value, generate_node_id
)


class TestAnnotationHypothesis:
    def test_to_dict(self):
        hyp = AnnotationHypothesis(
            cell_type="CL:0000084 (memory T cell)",
            cell_ontology_id="CL:0000084",
            reasoning="Marker profile matches",
            expected_markers=["CD3D", "CD4", "CCR7"],
            confidence_estimate=0.75
        )
        d = hyp.to_dict()
        assert d["cell_type"] == "CL:0000084 (memory T cell)"
        assert len(d["expected_markers"]) == 3


class TestMCTSNode:
    def test_initial_stats(self):
        hyp = AnnotationHypothesis("T cell", None, "test", [])
        node = MCTSNode("n1", "cluster_001", hyp)
        assert node.visits == 0
        assert node.avg_reward == 0.0
    
    def test_update_stats(self):
        hyp = AnnotationHypothesis("T cell", None, "test", [])
        node = MCTSNode("n1", "cluster_001", hyp)
        node.update(0.8, [], {"weighted_score": 0.8})
        assert node.visits == 1
        assert node.avg_reward == 0.8
        assert node.confidence == 0.8
    
    def test_ucb_unvisited(self):
        hyp = AnnotationHypothesis("T cell", None, "test", [])
        parent = MCTSNode("root", "c1", hyp, visits=10)
        child = MCTSNode("c1", "c1", hyp, parent=parent)
        assert child.ucb_score() == float('inf')
    
    def test_ucb_exploitation_exploration(self):
        hyp = AnnotationHypothesis("T cell", None, "test", [])
        parent = MCTSNode("root", "c1", hyp, visits=100)
        
        # High reward, many visits = high exploitation, low exploration
        child1 = MCTSNode("c1", "c1", hyp, parent=parent, visits=50, total_reward=40.0)
        # Low reward, few visits = low exploitation, high exploration
        child2 = MCTSNode("c2", "c1", hyp, parent=parent, visits=5, total_reward=2.0)
        
        ucb1 = child1.ucb_score(1.414, 100)
        ucb2 = child2.ucb_score(1.414, 100)
        
        # child1 should still win due to much higher avg reward
        assert ucb1 > ucb2
    
    def test_backpropagate(self):
        hyp = AnnotationHypothesis("T cell", None, "test", [])
        root = MCTSNode("root", "c1", hyp)
        child = MCTSNode("c1", "c1", hyp)
        root.add_child(child)
        
        backpropagate(child, 0.9)
        
        assert child.visits == 1
        assert child.total_reward == 0.9
        assert root.visits == 1  # Backpropagated
        assert root.total_reward == 0.9


class TestSelectionAndExtraction:
    def test_select_node_prefers_unvisited(self):
        hyp = AnnotationHypothesis("T", None, "r", [])
        root = MCTSNode("root", "c1", hyp, visits=10)
        
        visited = MCTSNode("v1", "c1", hyp, parent=root, visits=5, total_reward=4.0)
        unvisited = MCTSNode("u1", "c1", hyp, parent=root)
        root.children = [visited, unvisited]
        
        selected = select_node(root)
        assert selected == unvisited
    
    def test_extract_best(self):
        hyp = AnnotationHypothesis("T", None, "r", [])
        root = MCTSNode("root", "c1", hyp)
        
        # Add children with different rewards
        for i, reward in enumerate([0.3, 0.9, 0.5]):
            child = MCTSNode(f"c{i}", "c1", hyp, parent=root)
            child.update(reward, [], {"weighted_score": reward})
            root.add_child(child)
        
        best = extract_best(root)
        assert best.avg_reward == 0.9


class TestUtilityFunctions:
    def test_compute_annotation_value(self):
        evaluation = {
            "dimension_scores": {
                "marker_match": 0.9,
                "ontology_consistency": 0.8,
                "evidence_consistency": 0.7,
                "llm_confidence": 0.6
            }
        }
        weights = {
            "marker_match": 0.4,
            "ontology_consistency": 0.3,
            "evidence_consistency": 0.2,
            "llm_confidence": 0.1
        }
        score = compute_annotation_value(evaluation, weights)
        expected = 0.4*0.9 + 0.3*0.8 + 0.2*0.7 + 0.1*0.6
        assert abs(score - expected) < 0.001
    
    def test_generate_node_id_deterministic(self):
        id1 = generate_node_id("cluster_001", ["root", "T_cell"])
        id2 = generate_node_id("cluster_001", ["root", "T_cell"])
        assert id1 == id2
    
    def test_generate_node_id_unique_paths(self):
        id1 = generate_node_id("cluster_001", ["root", "T_cell"])
        id2 = generate_node_id("cluster_001", ["root", "B_cell"])
        assert id1 != id2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

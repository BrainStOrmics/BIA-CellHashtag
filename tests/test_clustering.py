"""Integration tests for clustering subgraph. Run with: pytest tests/test_clustering.py -v"""

import pickle
from pathlib import Path

import pytest
import scanpy as sc

from cellhashtag.graphs.clustering import (
    build_clustering_subgraph,
    _run_clustering,
    _assess_quality,
    _route_after_assessment,
)


@pytest.fixture
def mock_state(tmp_path):
    adata = sc.read_h5ad("data/example.h5ad")
    pkl_path = str(tmp_path / "test.pkl")
    with open(pkl_path, "wb") as f:
        pickle.dump(adata, f)
    return {
        "adata_path": pkl_path,
        "cluster_key": "leiden",
        "clustering_params": {"resolution": 0.5},
        "clustering_result": {},
        "clustering_quality": {},
        "clustering_history": [],
        "iteration_count": 0,
        "max_iterations": 2,
        "status": "starting",
        "errors": [],
    }


class TestClusteringSubgraph:
    def test_build_graph(self):
        assert build_clustering_subgraph() is not None

    def test_clustering_runs(self, mock_state):
        result = _run_clustering(mock_state)
        assert result["clustering_result"]["n_clusters"] > 0
        assert result["iteration_count"] == 1

    def test_quality_assessment(self, mock_state):
        mock_state.update(_run_clustering(mock_state))
        result = _assess_quality(mock_state)
        assert "clustering_quality" in result

    def test_routing_done_on_max_iterations(self, mock_state):
        mock_state["iteration_count"] = 10
        mock_state["max_iterations"] = 3
        assert _route_after_assessment(mock_state) == "done"

    def test_routing_retry_on_low_quality(self, mock_state):
        mock_state["iteration_count"] = 1
        mock_state["max_iterations"] = 5
        mock_state["clustering_quality"] = {"silhouette": {"recommendation": "adjust"}}
        assert _route_after_assessment(mock_state) == "retry"

    def test_routing_done_on_good_quality(self, mock_state):
        mock_state["iteration_count"] = 1
        mock_state["max_iterations"] = 5
        mock_state["clustering_quality"] = {
            "silhouette": {"recommendation": "ok"},
            "modularity": {"recommendation": "ok"},
        }
        assert _route_after_assessment(mock_state) == "done"

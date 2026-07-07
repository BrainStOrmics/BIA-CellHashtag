"""Integration tests for orchestrator graph. Run with: pytest tests/test_orchestrator.py -v"""

import pickle
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import scanpy as sc

from cellhashtag.graphs.orchestrator import (
    build_orchestrator_graph,
    node_perception,
    node_clustering,
    node_annotation,
    node_output,
    _parse_agent_result,
)


@pytest.fixture
def mock_state(tmp_path):
    adata = sc.read_h5ad("data/example.h5ad")
    pkl_path = str(tmp_path / "test.pkl")
    with open(pkl_path, "wb") as f:
        pickle.dump(adata, f)
    return {
        "adata_path": str(tmp_path / "test.h5ad"),
        "input_format": "h5ad",
        "omics_type": "auto",
        "cluster_key": "leiden",
        "markers_per_cluster": {},
        "annotation_results": [],
        "output_dir": str(tmp_path / "output"),
        "status": "starting",
        "errors": [],
    }


class TestOrchestratorGraph:
    def test_build_graph(self):
        graph = build_orchestrator_graph()
        assert graph is not None

    def test_graph_structure(self):
        graph = build_orchestrator_graph()
        nodes = graph.nodes
        assert "perception" in nodes
        assert "clustering" in nodes
        assert "annotation" in nodes
        assert "output" in nodes


class TestPerceptionNode:
    @patch("scanpy.pp.filter_cells")
    @patch("scanpy.pp.neighbors")
    @patch("scanpy.pp.pca")
    @patch("scanpy.pp.highly_variable_genes")
    def test_perception_runs(self, mock_hvg, mock_pca, mock_neighbors, mock_filter, mock_state, tmp_path):
        # Copy h5ad to expected location
        import shutil
        shutil.copy("data/example.h5ad", mock_state["adata_path"])

        result = node_perception(mock_state)
        assert "adata_path" in result
        assert result["adata_path"].endswith(".pkl")
        assert "omics_type" in result
        assert "markers_per_cluster" in result
        # Mock data has only 150 genes, so markers may be empty after filtering


class TestClusteringNode:
    @patch("cellhashtag.config.config.load_config")
    def test_clustering_runs(self, mock_load_config, mock_state, tmp_path):
        # Setup mock config
        from cellhashtag.config import CellHashtagConfig
        mock_load_config.return_value = CellHashtagConfig()

        # Create pickle file
        adata = sc.read_h5ad("data/example.h5ad")
        with open(mock_state["adata_path"].replace(".h5ad", ".pkl"), "wb") as f:
            pickle.dump(adata, f)
        mock_state["adata_path"] = mock_state["adata_path"].replace(".h5ad", ".pkl")

        result = node_clustering(mock_state)
        assert "status" in result
        assert "clustering_done" in result["status"]


class TestAnnotationNode:
    @patch("cellhashtag.config.setup_llm")
    @patch("deepagents.create_deep_agent")
    @patch("cellhashtag.config.config.load_config")
    def test_annotation_runs(
        self, mock_load_config, mock_create_agent, mock_setup_llm, mock_state, tmp_path
    ):
        # Setup mocks
        from cellhashtag.config import CellHashtagConfig
        mock_load_config.return_value = CellHashtagConfig()
        mock_setup_llm.return_value = None

        # Mock DeepAgents response
        mock_agent = Mock()
        mock_agent.invoke.return_value = {
            "messages": [
                Mock(content='{"cell_type": "T cell", "confidence": 0.9, "reasoning": "High CD3D"}')
            ]
        }
        mock_create_agent.return_value = mock_agent

        # Create pickle file with clusters
        adata = sc.read_h5ad("data/example.h5ad")
        pkl_path = mock_state["adata_path"].replace(".h5ad", ".pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(adata, f)
        mock_state["adata_path"] = pkl_path
        mock_state["markers_per_cluster"] = {"0": ["CD3D", "CD3E"], "1": ["CD19", "MS4A1"]}

        result = node_annotation(mock_state)
        assert "annotation_results" in result
        assert len(result["annotation_results"]) > 0
        assert result["status"] == "annotation_done"


class TestOutputNode:
    def test_output_saves_files(self, mock_state, tmp_path):
        # Create pickle file with annotations
        adata = sc.read_h5ad("data/example.h5ad")
        pkl_path = tmp_path / "test.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(adata, f)
        mock_state["adata_path"] = str(pkl_path)
        mock_state["annotation_results"] = [
            {"cluster": "0", "cell_type": "T cell", "confidence": 0.9, "reasoning": "test"},
            {"cluster": "1", "cell_type": "B cell", "confidence": 0.85, "reasoning": "test"},
        ]

        result = node_output(mock_state)
        assert result["status"] == "done"

        # Check files created
        output_dir = Path(mock_state["output_dir"])
        assert (output_dir / "annotated.h5ad").exists()
        assert (output_dir / "annotation_table.csv").exists()


class TestParseAgentResult:
    def test_parse_valid_json(self):
        result = {
            "messages": [
                Mock(content='{"cell_type": "T cell", "confidence": 0.95, "reasoning": "CD3D high"}')
            ]
        }
        parsed = _parse_agent_result(result, "cluster_0")
        assert parsed["cell_type"] == "T cell"
        assert parsed["confidence"] == 0.95

    def test_parse_no_json(self):
        result = {"messages": [Mock(content="No JSON here")]}
        parsed = _parse_agent_result(result, "cluster_0")
        assert parsed["cell_type"] == "Unknown"
        assert parsed["confidence"] == 0.3

    def test_parse_invalid_json(self):
        result = {"messages": [Mock(content='{"cell_type": "T cell"')]}  # Malformed
        parsed = _parse_agent_result(result, "cluster_0")
        assert parsed["cell_type"] == "Unknown"

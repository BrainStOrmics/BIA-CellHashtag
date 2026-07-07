"""CellHashtag v3.0 — Entry point."""

from pathlib import Path
from typing import Any, Optional

from langgraph.checkpoint.memory import MemorySaver

from .config import load_config, setup_llm, CellHashtagConfig, PROFILES
from .graphs import build_orchestrator_graph


class CellHashtagAgent:
    def __init__(
        self,
        config_path: Optional[str] = None,
        *,
        profile: str = "default",
        **config_overrides,
    ):
        """Initialize CellHashtag Agent.

        Args:
            config_path: Path to YAML config. None → default.
            profile: Preset mode — "fast" (quick/cheap), "default", "deep" (thorough).
            **config_overrides: Runtime config overrides using __ for nesting,
                e.g. llm__model="gpt-4o", annotation__max_anno_iter=3.
        """
        self.config = load_config(config_path, profile=profile, **config_overrides)
        self.graph = build_orchestrator_graph().compile(checkpointer=MemorySaver())

    def run(
        self,
        input_path: str,
        output_dir: str = "output",
        input_format: str = "h5ad",
        omics_type: str = "auto",
        cluster_key: str = "leiden",
    ) -> dict:
        state = {
            "adata_path": input_path,
            "input_format": input_format,
            "omics_type": omics_type,
            "cluster_key": cluster_key,
            "markers_per_cluster": [],
            "annotation_results": [],
            "output_dir": output_dir,
            "status": "starting",
            "errors": [],
        }
        result = self.graph.invoke(state)
        return {
            "output_dir": output_dir,
            "n_clusters": len(result.get("annotation_results", [])),
            "status": result.get("status", "unknown"),
            "errors": result.get("errors", []),
        }

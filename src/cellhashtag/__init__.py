"""
CellHashtag v2.0 — LLM-Based Single-Cell Annotation Agent

Automated cell type annotation for scRNA-seq and spatial RNA-seq data.
Features:
- Built-in clustering with quality assessment (Phiclust, Silhouette, Modularity)
- Multi-format input (h5ad, rds, 10x mtx)
- LLM self-criticism annotation loop
- CellWiki knowledge integration
- Dual-layer annotation (Cell Type + Sub-annotation)
- HITL interaction support
"""

from .agent import CellHashtagAgent
from .state import CellHashtagState
from .config.config import LLMConfig, AgentConfig, load_config, setup_llm

__all__ = [
    "CellHashtagAgent",
    "CellHashtagState",
    "LLMConfig",
    "AgentConfig",
    "load_config",
    "setup_llm",
]

__version__ = "2.0.0"

"""CellHashtag v3.0 — LLM-Based Single-Cell Annotation Agent."""

from .agent import CellHashtagAgent
from .config.config import (
    CellHashtagConfig,
    LLMConfig,
    ClusteringConfig,
    AnnotationConfig,
    LATSConfig,
    load_config,
    setup_llm,
    PROFILES,
)

__all__ = [
    "CellHashtagAgent",
    "CellHashtagConfig",
    "LLMConfig",
    "ClusteringConfig",
    "AnnotationConfig",
    "LATSConfig",
    "load_config",
    "setup_llm",
    "PROFILES",
]

__version__ = "3.0.0.dev"

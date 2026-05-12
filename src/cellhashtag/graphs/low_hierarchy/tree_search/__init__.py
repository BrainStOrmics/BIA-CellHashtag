"""
LATS Tree Search Subgraph for CellHashtag

Exports:
    create_graph(): Build the LATS search LangGraph
    LATSState: TypedDict for graph state
"""
from typing import TypedDict, List, Optional, Dict, Any
from .graph import create_graph

__all__ = ["create_graph", "LATSState"]


class LATSState(TypedDict):
    """State for LATS tree search graph"""
    # Context
    cluster_id: str
    cluster_markers: List[str]
    tissue_source: str
    experimental_condition: Optional[str]
    
    # Search state
    current_hypotheses: List[Dict]  # List of AnnotationHypothesis dicts
    evidence_cache: Dict[str, List[Dict]]  # marker -> evidence list
    search_tree_root: Optional[Dict]  # Lightweight MCTS root summary
    
    # Iteration control
    iteration: int
    max_iterations: int
    confidence_threshold: float
    
    # Results
    best_annotation: Optional[Dict]
    best_confidence: float
    search_complete: bool
    search_status: str  # "continue" | "confident" | "exhausted"
    
    # Config references (paths, not objects)
    config_path: str
    prompts_dir: str

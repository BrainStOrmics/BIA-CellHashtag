"""
Low Hierarchy Annotation Graph for CellHashtag

Implements the swarm-mode annotation with LATS tree search as default strategy.

Exports:
    create_graph(): Build the low hierarchy LangGraph
    LowHierarchyState: TypedDict for graph state
"""
from typing import TypedDict, List, Optional, Dict, Any
from .graph import create_graph

# Re-export tree_search components for convenience
from .tree_search import LATSState as TreeSearchState
from .tree_search import create_graph as create_tree_search_graph

__all__ = ["create_graph", "LowHierarchyState", "create_tree_search_graph", "TreeSearchState"]


class LowHierarchyState(TypedDict):
    """
    State for low hierarchy annotation graph.
    
    Combines cluster context with LATS search state.
    """
    # Cluster context (from high hierarchy)
    cluster_id: str
    cluster_markers: List[str]
    cluster_expression_summary: Dict[str, float]  # gene -> avg expression
    tissue_source: str
    experimental_condition: Optional[str]
    
    # Annotation constraints (from Harness node)
    allowed_cell_types: Optional[List[str]]  # Cell Ontology IDs to consider
    excluded_markers: List[str]  # Markers that should NOT be highly expressed
    priority_pathways: List[str]  # Biological pathways to prioritize
    
    # LATS search state (delegated to tree_search subgraph)
    lats_config_path: str  # Path to config.yaml
    lats_prompts_dir: str  # Path to prompts/low_hierarchy
    search_result: Optional[Dict]  # Output from tree_search subgraph
    
    # Iteration and control
    annotation_iteration: int
    max_annotation_iter: int
    min_confidence: float  # Minimum confidence to accept annotation
    
    # Output
    annotation_result: Optional[Dict]  # Final {cell_type, confidence, evidence}
    annotation_complete: bool
    fallback_to_manual: bool  # If search exhausted without confident result

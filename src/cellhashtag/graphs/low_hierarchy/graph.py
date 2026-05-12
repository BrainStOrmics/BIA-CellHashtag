"""
Low Hierarchy Annotation Graph

Default implementation using LATS tree search for annotation.
Workflow: Harness(constraints) → LATS_Search → Memory(record) → Return
"""
from langgraph.graph import StateGraph, END
from typing import Dict, Any, cast

from . import LowHierarchyState
from .tree_search import create_tree_search_graph, TreeSearchState
from .nodes import apply_constraints, record_to_memory, finalize_annotation


def create_graph(global_config: Dict[str, Any]) -> StateGraph:
    """
    Build the low hierarchy annotation LangGraph with LATS as default search.
    
    Args:
        global_config: Full config dict from config.yaml
        
    Returns:
        Compiled StateGraph ready for execution
    """
    # Initialize main workflow
    workflow = StateGraph(LowHierarchyState)
    
    # Add nodes
    workflow.add_node("harness", apply_constraints)
    workflow.add_node("lats_search", run_lats_subgraph)
    workflow.add_node("memory", record_to_memory)
    workflow.add_node("finalize", finalize_annotation)
    
    # Define edges
    workflow.set_entry_point("harness")
    workflow.add_edge("harness", "lats_search")
    workflow.add_edge("lats_search", "memory")
    workflow.add_edge("memory", "finalize")
    workflow.add_edge("finalize", END)
    
    return workflow


def run_lats_subgraph(state: LowHierarchyState) -> Dict:
    """
    Execute LATS tree search as a subgraph.
    
    Maps LowHierarchyState → TreeSearchState → results → LowHierarchyState
    """
    from langgraph.graph import StateGraph
    
    # Prepare LATS subgraph state
    lats_state: TreeSearchState = {
        # Context mapping
        "cluster_id": state["cluster_id"],
        "cluster_markers": state["cluster_markers"],
        "tissue_source": state["tissue_source"],
        "experimental_condition": state.get("experimental_condition"),
        
        # Initialize search state
        "current_hypotheses": [],
        "evidence_cache": {},
        "search_tree_root": None,
        "iteration": 0,
        
        # Config from global config
        "max_iterations": global_config.get("lats", {}).get("search_params", {}).get("n_iterations", 10),
        "confidence_threshold": global_config.get("lats", {}).get("search_params", {}).get("confidence_threshold", 0.7),
        "best_confidence": 0.0,
        "search_complete": False,
        "search_status": "continue",
        
        # File paths for loading prompts/config
        "config_path": state["lats_config_path"],
        "prompts_dir": state["lats_prompts_dir"],
    }
    
    # Build and compile LATS subgraph
    lats_graph = create_tree_search_graph(global_config)
    compiled_lats = lats_graph.compile()
    
    # Execute search
    result = compiled_lats.invoke(lats_state)
    
    # Map results back to low hierarchy state
    search_result = result.get("best_annotation") if result else None
    
    return {
        "search_result": search_result,
        # Pass through confidence for finalization decision
        "_lats_confidence": result.get("best_confidence", 0.0) if result else 0.0,
        "_lats_iterations": result.get("total_iterations", 0) if result else 0,
    }


def route_by_confidence(state: LowHierarchyState) -> str:
    """
    Route based on annotation confidence.
    
    - High confidence: finalize and return
    - Low confidence + iterations left: could retry (future extension)
    - Exhausted: mark for manual review
    """
    confidence = state.get("_lats_confidence", 0.0)
    min_conf = state.get("min_confidence", 0.7)
    
    if confidence >= min_conf:
        return "accept"
    elif state.get("annotation_iteration", 0) < state.get("max_annotation_iter", 3):
        return "retry"  # Future: implement retry logic
    else:
        return "fallback"

"""
LATS Tree Search Graph Construction

Builds a LangGraph for MCTS-based annotation search.
"""
from langgraph.graph import StateGraph, END
from typing import Dict, Any

from . import LATSState
from .nodes import (
    initialize_search,
    select_and_expand,
    validate_hypotheses,
    evaluate_and_backpropagate,
    check_termination,
    extract_result
)


def create_graph(config: Dict[str, Any]) -> StateGraph:
    """
    Build the LATS tree search LangGraph.
    
    Args:
        config: Configuration dict with lats.search_params, etc.
    
    Returns:
        Compiled StateGraph ready for execution
    """
    # Initialize graph with state schema
    workflow = StateGraph(LATSState)
    
    # Add nodes
    workflow.add_node("init", initialize_search)
    workflow.add_node("select_expand", select_and_expand)
    workflow.add_node("validate", validate_hypotheses)
    workflow.add_node("evaluate", evaluate_and_backpropagate)
    workflow.add_node("check_done", check_termination)
    workflow.add_node("extract", extract_result)
    
    # Define edges
    workflow.set_entry_point("init")
    workflow.add_edge("init", "select_expand")
    workflow.add_edge("select_expand", "validate")
    workflow.add_edge("validate", "evaluate")
    workflow.add_edge("evaluate", "check_done")
    
    # Conditional routing from check_done
    workflow.add_conditional_edges(
        "check_done",
        route_after_evaluation,
        {
            "continue": "select_expand",
            "confident": "extract",
            "exhausted": "extract",
        }
    )
    
    workflow.add_edge("extract", END)
    
    return workflow


def route_after_evaluation(state: LATSState) -> str:
    """Route based on search status after evaluation"""
    status = state.get("search_status", "continue")
    
    if status == "confident":
        return "confident"
    elif status == "exhausted":
        return "exhausted"
    else:
        return "continue"

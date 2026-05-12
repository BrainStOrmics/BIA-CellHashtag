"""
Low Hierarchy Node Implementations

Nodes for constraint application, memory recording, and result finalization.
LATS search is delegated to tree_search subgraph.
"""
from typing import Dict, List, Any, Optional
import json
import os

from ..utils.utilities import load_config
from ..utils.cellwiki_client import query_cell_type_constraints


def apply_constraints(state: Dict) -> Dict:
    """
    Harness node: Apply annotation constraints before search.
    
    Filters allowed cell types based on tissue context and excludes
    incompatible marker patterns.
    
    Input state keys:
        tissue_source, cluster_markers, allowed_cell_types, excluded_markers
    
    Returns:
        Constrained search parameters for LATS
    """
    tissue = state.get("tissue_source", "")
    cluster_markers = set(state.get("cluster_markers", []))
    
    # Default: allow all cell types if not specified
    allowed = state.get("allowed_cell_types")
    if allowed is None:
        # Query CellWiki for tissue-appropriate cell types
        tissue_types = query_cell_type_constraints(tissue)
        allowed = tissue_types.get("compatible_types", [])
    
    # Filter by excluded markers
    if state.get("excluded_markers"):
        excluded = set(state["excluded_markers"])
        # In production: query which cell types express excluded markers
        # For now, just log the constraint
        print(f"Constraint: exclude types expressing {excluded & cluster_markers}")
    
    # Add priority pathways as search hints
    priority = state.get("priority_pathways", [])
    search_hints = {
        "focus_markers": [m for m in cluster_markers if any(p.lower() in m.lower() for p in priority)],
        "tissue_context": tissue,
        "allowed_ontology_ids": allowed[:20] if allowed else None  # Limit for prompt context
    }
    
    return {
        "lats_search_hints": search_hints,
        "_constraints_applied": True
    }


def record_to_memory(state: Dict) -> Dict:
    """
    Memory node: Record search results for future reference.
    
    Caches annotation results to avoid redundant searches for similar clusters.
    
    Returns:
        Updated state with memory entry confirmation
    """
    search_result = state.get("search_result")
    if not search_result:
        return {"_memory_recorded": False, "memory_key": None}
    
    # Create cache key from cluster signature
    cluster_id = state["cluster_id"]
    markers_sig = tuple(sorted(state["cluster_markers"][:10]))
    tissue = state["tissue_source"]
    memory_key = f"{tissue}:{hash(markers_sig) & 0xFFFF:04x}"
    
    # In production: write to persistent cache (Redis/SQLite)
    # For now: store in state for this session
    memory_entry = {
        "cluster_id": cluster_id,
        "annotation": search_result.get("cell_type"),
        "confidence": search_result.get("confidence"),
        "evidence_summary": [e.get("content", "")[:100] for e in search_result.get("evidence", [])[:3]],
        "timestamp": "__session__"  # Would be actual timestamp in production
    }
    
    # Append to session memory list
    existing_memory = state.get("_session_memory", [])
    existing_memory.append((memory_key, memory_entry))
    
    return {
        "_memory_recorded": True,
        "memory_key": memory_key,
        "_session_memory": existing_memory
    }


def finalize_annotation(state: Dict) -> Dict:
    """
    Finalize node: Format and validate final annotation result.
    
    Applies confidence threshold and formats output for high hierarchy.
    
    Returns:
        annotation_result with standardized schema
    """
    search_result = state.get("search_result")
    confidence = state.get("_lats_confidence", 0.0)
    min_conf = state.get("min_confidence", 0.7)
    
    # Handle no result case
    if not search_result:
        return {
            "annotation_result": None,
            "annotation_complete": True,
            "fallback_to_manual": True,
            "fallback_reason": "LATS search returned no valid annotation"
        }
    
    # Check confidence threshold
    if confidence < min_conf:
        return {
            "annotation_result": {
                "cell_type": search_result.get("cell_type"),
                "confidence": confidence,
                "status": "low_confidence",
                "evidence": search_result.get("evidence", []),
                "recommendation": "Manual review recommended"
            },
            "annotation_complete": True,
            "fallback_to_manual": True,
            "fallback_reason": f"Confidence {confidence:.2f} < threshold {min_conf}"
        }
    
    # Format successful result
    result = {
        "cell_type": search_result["cell_type"],
        "cell_ontology_id": extract_ontology_id(search_result["cell_type"]),
        "confidence": confidence,
        "evidence": search_result.get("evidence", [])[:5],  # Top 5 evidence items
        "reasoning": search_result.get("reasoning", []),
        "search_metadata": {
            "iterations": state.get("_lats_iterations", 0),
            "memory_key": state.get("memory_key")
        }
    }
    
    return {
        "annotation_result": result,
        "annotation_complete": True,
        "fallback_to_manual": False
    }


def extract_ontology_id(cell_type_str: str) -> Optional[str]:
    """Extract Cell Ontology ID from formatted string like 'CL:0000084 (T cell)'"""
    if not cell_type_str:
        return None
    # Pattern: CL:XXXXXXX or http://purl.obolibrary.org/obo/CL_XXXXXXX
    if cell_type_str.startswith("CL:"):
        return cell_type_str.split()[0]
    if "obo/CL_" in cell_type_str:
        return "CL:" + cell_type_str.split("CL_")[-1].split()[0].replace(")", "")
    return None

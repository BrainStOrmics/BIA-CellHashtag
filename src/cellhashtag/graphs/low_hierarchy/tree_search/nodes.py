"""
LATS Tree Search Node Implementations

Each function is a LangGraph node that processes state and returns updates.
"""
import os
import json
import uuid
from typing import Dict, List, Any, Optional

from ...skills.lats_search.mcts_core import (
    MCTSNode, AnnotationHypothesis, EvidenceItem,
    select_node, backpropagate, extract_best,
    compute_annotation_value, generate_node_id
)
from ..utils.llm_client import call_llm_with_prompt
from ..utils.cellwiki_client import query_cellwiki
from ..utils.utilities import load_config, load_prompt_template


def initialize_search(state: Dict) -> Dict:
    """
    Initialize MCTS search tree from cluster context.
    
    Input state keys:
        cluster_id, cluster_markers, tissue_source, config_path, prompts_dir
    
    Returns:
        Updated state with initialized search_tree_root and iteration=0
    """
    config = load_config(state["config_path"])
    lats_config = config.get("lats", {})
    
    # Create root hypothesis from initial context
    root_hypothesis = AnnotationHypothesis(
        cell_type="UNASSIGNED",
        cell_ontology_id=None,
        reasoning="Root node representing unannotated cluster",
        expected_markers=state["cluster_markers"][:5],  # Top 5 markers
        confidence_estimate=0.0
    )
    
    root = MCTSNode(
        node_id=generate_node_id(state["cluster_id"], ["root"]),
        cluster_id=state["cluster_id"],
        hypothesis=root_hypothesis,
        visits=1,  # Root is always "visited"
        total_reward=0.0
    )
    
    return {
        "search_tree_root": root.to_summary(),
        "iteration": 0,
        "max_iterations": lats_config.get("search_params", {}).get("n_iterations", 10),
        "confidence_threshold": lats_config.get("search_params", {}).get("confidence_threshold", 0.7),
        "best_confidence": 0.0,
        "search_complete": False,
        "search_status": "continue",
        "evidence_cache": {}
    }


def select_and_expand(state: Dict) -> Dict:
    """
    Select node via UCB and expand with LLM-generated hypotheses.
    
    Returns:
        New hypotheses to validate, updated tree summary
    """
    from ...skills.lats_search.mcts_core import MCTSNode
    
    # Reconstruct minimal tree from summary (in production, use persistent storage)
    root_summary = state["search_tree_root"]
    if not root_summary:
        return {"new_hypotheses": [], "expansion_failed": True}
    
    # For now, create a simple root for selection
    # In full implementation: deserialize tree from storage
    root = _reconstruct_root_from_summary(root_summary, state)
    
    # Select node for expansion
    lats_config = load_config(state["config_path"]).get("lats", {})
    exploration_weight = lats_config.get("search_params", {}).get("exploration_weight", 1.414)
    max_branches = lats_config.get("search_params", {}).get("max_branches", 3)
    
    selected = select_node(root, exploration_weight)
    
    # Skip expansion if already fully expanded or root with no progress
    if selected.is_fully_expanded(max_branches) or (selected == root and selected.visits > 1):
        return {"new_hypotheses": [], "node_selected": selected.to_summary()}
    
    # Generate new hypotheses via LLM
    prompt = load_prompt_template(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "prompts", "low_hierarchy", "lats_expansion.md"),
        cluster_id=state["cluster_id"],
        tissue_source=state["tissue_source"],
        experimental_condition=state.get("experimental_condition", "N/A"),
        top_markers_list=", ".join(state["cluster_markers"][:10]),
        current_candidates=[h.get("cell_type") for h in state.get("current_hypotheses", [])],
        evidence_summary=_summarize_evidence(state.get("evidence_cache", {}))
    )
    
    llm_config = lats_config.get("llm_config", {})
    response = call_llm_with_prompt(
        prompt,
        model=llm_config.get("model", "qwen3.5-plus"),
        temperature=llm_config.get("temperature", 0.2),
        response_format="json_object"
    )
    
    try:
        result = json.loads(response)
        hypotheses_data = result.get("hypotheses", [])
    except json.JSONDecodeError:
        hypotheses_data = []
    
    # Create child nodes for new hypotheses
    new_hypotheses = []
    for hyp_data in hypotheses_data[:max_branches]:
        hyp = AnnotationHypothesis(
            cell_type=hyp_data.get("cell_type", "UNKNOWN"),
            cell_ontology_id=hyp_data.get("cell_ontology_id"),
            reasoning=hyp_data.get("reasoning", ""),
            expected_markers=hyp_data.get("expected_markers", []),
            confidence_estimate=hyp_data.get("confidence_estimate", 0.5),
            distinguishing_features=hyp_data.get("distinguishing_features", "")
        )
        
        child = MCTSNode(
            node_id=generate_node_id(state["cluster_id"], [selected.node_id, hyp.cell_type]),
            cluster_id=state["cluster_id"],
            hypothesis=hyp,
            parent=selected
        )
        selected.add_child(child)
        new_hypotheses.append(hyp.to_dict())
    
    return {
        "new_hypotheses": new_hypotheses,
        "search_tree_root": root.to_summary(),
        "selected_node": selected.to_summary()
    }


def validate_hypotheses(state: Dict) -> Dict:
    """
    Parallel validation of new hypotheses against evidence sources.
    
    Returns:
        Validation results for each hypothesis
    """
    hypotheses = state.get("new_hypotheses", [])
    if not hypotheses:
        return {"validation_results": []}
    
    results = []
    evidence_cache = state.get("evidence_cache", {})
    
    for hyp in hypotheses:
        expected_markers = hyp.get("expected_markers", [])
        cell_type = hyp.get("cell_type")
        
        # Query CellWiki for marker validation
        cw_evidence = []
        for marker in expected_markers[:5]:  # Limit API calls
            cache_key = f"{cell_type}:{marker}"
            if cache_key in evidence_cache:
                cw_evidence.extend(evidence_cache[cache_key])
            else:
                marker_info = query_cellwiki(marker, cell_type)
                evidence = [
                    EvidenceItem(
                        source="cellwiki",
                        content=marker_info.get("description", ""),
                        confidence=marker_info.get("confidence", 0.7),
                        supports=marker_info.get("matches_cluster", True)
                    )
                ]
                cw_evidence.extend(evidence)
                evidence_cache[cache_key] = [e.to_dict() for e in evidence]
        
        # Compute marker overlap score
        cluster_markers = set(state.get("cluster_markers", []))
        expected_set = set(expected_markers)
        overlap = len(cluster_markers & expected_set) / max(len(expected_set), 1)
        
        results.append({
            "hypothesis": hyp,
            "marker_overlap": overlap,
            "evidence_items": [e.to_dict() for e in cw_evidence],
            "validation_summary": f"Overlap: {overlap:.2f}, Evidence count: {len(cw_evidence)}"
        })
    
    return {
        "validation_results": results,
        "evidence_cache": evidence_cache
    }


def evaluate_and_backpropagate(state: Dict) -> Dict:
    """
    Evaluate validation results and backpropagate rewards through tree.
    
    Returns:
        Updated best confidence, evaluation details
    """
    from ...skills.lats_search.mcts_core import MCTSNode, compute_annotation_value
    
    results = state.get("validation_results", [])
    if not results:
        return {"evaluations": []}
    
    config = load_config(state["config_path"])
    weights = config.get("lats", {}).get("value_function", {}).get("weights", {
        "marker_match": 0.4,
        "ontology_consistency": 0.3,
        "evidence_consistency": 0.2,
        "llm_confidence": 0.1
    })
    
    evaluations = []
    lats_config = config.get("lats", {})
    llm_config = lats_config.get("llm_config", {})
    
    for res in results:
        # Prepare evaluation prompt
        hyp = res["hypothesis"]
        prompt = load_prompt_template(
            os.path.join(os.path.dirname(__file__), "..", "..", "..", "prompts", "low_hierarchy", "lats_evaluation.md"),
            candidate_cell_type=hyp.get("cell_type"),
            tissue_source=state["tissue_source"],
            cluster_top_markers=", ".join(state["cluster_markers"][:10]),
            cluster_quality_score=0.85,  # Would come from clustering_quality skill
            supporting_evidence_list="\n".join([
                f"- [{e['source']}] {e['content']}" 
                for e in res["evidence_items"] if e.get("supports")
            ]),
            conflicting_evidence_list="\n".join([
                f"- [{e['source']}] {e['content']}" 
                for e in res["evidence_items"] if not e.get("supports")
            ]),
            expected_markers=", ".join(hyp.get("expected_markers", [])),
            observed_overlap=f"{res['marker_overlap']:.2%}"
        )
        
        # Get LLM evaluation
        response = call_llm_with_prompt(
            prompt,
            model=llm_config.get("model", "qwen3.5-plus"),
            temperature=0.1,  # Very low for consistent evaluation
            response_format="json_object"
        )
        
        try:
            evaluation = json.loads(response)
        except json.JSONDecodeError:
            evaluation = {"weighted_score": 0.5, "recommendation": "need_more_evidence"}
        
        # Compute final value
        weighted_score = compute_annotation_value(evaluation, weights)
        
        evaluations.append({
            "hypothesis_id": hyp.get("cell_type"),
            "evaluation": evaluation,
            "weighted_score": weighted_score,
            "recommendation": evaluation.get("recommendation", "need_more_evidence")
        })
        
        # Backpropagate through tree (simplified: update root stats)
        # In full implementation: traverse from leaf to root
        root_summary = state.get("search_tree_root", {})
        if root_summary:
            current_best = max(
                weighted_score, 
                root_summary.get("avg_reward", 0.0)
            )
            root_summary["avg_reward"] = current_best
            root_summary["confidence"] = max(
                root_summary.get("confidence", 0.0),
                weighted_score
            )
    
    # Update best confidence
    best_score = max((e["weighted_score"] for e in evaluations), default=0.0)
    current_best = state.get("best_confidence", 0.0)
    
    return {
        "evaluations": evaluations,
        "best_confidence": max(current_best, best_score),
        "search_tree_root": state.get("search_tree_root")
    }


def check_termination(state: Dict) -> Dict:
    """
    Check if search should terminate based on confidence or iteration limits.
    
    Returns:
        search_status: "confident" | "exhausted" | "continue"
    """
    iteration = state.get("iteration", 0)
    max_iter = state.get("max_iterations", 10)
    confidence = state.get("best_confidence", 0.0)
    threshold = state.get("confidence_threshold", 0.7)
    
    # Check early stop threshold
    early_stop = load_config(state["config_path"]).get("lats", {}).get(
        "search_params", {}
    ).get("early_stop_threshold", 0.95)
    
    if confidence >= early_stop:
        status = "confident"
    elif confidence >= threshold or iteration >= max_iter:
        status = "exhausted" if iteration >= max_iter else "confident"
    else:
        status = "continue"
    
    return {
        "search_status": status,
        "iteration": iteration + 1,
        "search_complete": status != "continue"
    }


def extract_result(state: Dict) -> Dict:
    """
    Extract final annotation result from search tree.
    
    Returns:
        best_annotation with confidence and evidence
    """
    evaluations = state.get("evaluations", [])
    if not evaluations:
        return {
            "best_annotation": None,
            "best_confidence": 0.0,
            "search_complete": True
        }
    
    # Find best evaluation
    best_eval = max(evaluations, key=lambda e: e["weighted_score"])
    best_hyp = best_eval["hypothesis_id"]
    
    # Gather supporting evidence
    validation_results = state.get("validation_results", [])
    supporting_evidence = []
    for vr in validation_results:
        if vr["hypothesis"].get("cell_type") == best_hyp:
            supporting_evidence = vr["evidence_items"]
            break
    
    return {
        "best_annotation": {
            "cell_type": best_hyp,
            "confidence": best_eval["weighted_score"],
            "reasoning": best_eval["evaluation"].get("key_supporting_points", []),
            "evidence": supporting_evidence[:5],  # Top 5 evidence items
            "next_steps": best_eval["evaluation"].get("next_validation_suggestion")
        },
        "best_confidence": best_eval["weighted_score"],
        "search_complete": True,
        "total_iterations": state.get("iteration", 0)
    }


# Helper functions

def _reconstruct_root_from_summary(summary: Dict, state: Dict) -> MCTSNode:
    """Reconstruct minimal MCTSNode from lightweight summary"""
    from ...skills.lats_search.mcts_core import MCTSNode, AnnotationHypothesis
    
    hyp = AnnotationHypothesis(
        cell_type=summary.get("cell_type", "UNASSIGNED"),
        cell_ontology_id=None,
        reasoning="Reconstructed from summary",
        expected_markers=state["cluster_markers"][:5]
    )
    
    return MCTSNode(
        node_id=summary.get("node_id", "root"),
        cluster_id=state["cluster_id"],
        hypothesis=hyp,
        visits=summary.get("visits", 1),
        total_reward=summary.get("avg_reward", 0.0) * summary.get("visits", 1)
    )


def _summarize_evidence(evidence_cache: Dict) -> str:
    """Create brief summary of cached evidence for prompt context"""
    if not evidence_cache:
        return "No prior evidence collected."
    
    items = []
    for key, evidence_list in list(evidence_cache.items())[:3]:
        for e in evidence_list[:2]:
            support_marker = "✓" if e.get("supports") else "✗"
            items.append(f"{support_marker}[{e.get('source')}]: {e.get('content')[:50]}...")
    
    return "\n".join(items) if items else "Evidence pending collection."

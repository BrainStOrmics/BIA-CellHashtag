"""CellHashtag core algorithms (DAG-MCTS era)."""

from .dag_mcts import (
    DAGNode,
    TranspositionTable,
    run_dag_mcts,
    make_state_hash,
    select_node,
    expand,
    evaluate,
    backpropagate,
    extract_best,
    priority,
    fpu,
    soft_bellman_v,
    shared_discount,
    temperature_at_depth,
)
from .reflection import ReflectionSchema, parse_reflection
from .reward import (
    ClusterStats,
    RewardBreakdown,
    SlidingWindowNormalizer,
    compute_reward,
    final_reward,
    potential,
)
from .cold_start import ColdStartResult, cold_start
from .tree_viz import TreeVisualizer
from .lats_loop import run_lats_search, DEFAULT_WEIGHTS, DEFAULT_SEARCH_PARAMS
from .cache import EvidenceCache

__all__ = [
    "DAGNode",
    "TranspositionTable",
    "run_dag_mcts",
    "make_state_hash",
    "select_node",
    "expand",
    "evaluate",
    "backpropagate",
    "extract_best",
    "priority",
    "fpu",
    "soft_bellman_v",
    "shared_discount",
    "temperature_at_depth",
    "ReflectionSchema",
    "parse_reflection",
    "ClusterStats",
    "RewardBreakdown",
    "SlidingWindowNormalizer",
    "compute_reward",
    "final_reward",
    "potential",
    "ColdStartResult",
    "cold_start",
    "TreeVisualizer",
    "run_lats_search",
    "DEFAULT_WEIGHTS",
    "DEFAULT_SEARCH_PARAMS",
    "EvidenceCache",
]

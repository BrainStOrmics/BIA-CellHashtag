"""
MCTS Core Utilities for LATS-based Cell Annotation
"""
import math
import json
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
from enum import Enum


class SearchStatus(Enum):
    CONTINUE = "continue"
    CONFIDENT = "confident"
    EXHAUSTED = "exhausted"


@dataclass
class AnnotationHypothesis:
    """Represents a cell type annotation hypothesis"""
    cell_type: str
    cell_ontology_id: Optional[str]
    reasoning: str
    expected_markers: List[str]
    confidence_estimate: float = 0.5
    distinguishing_features: str = ""
    
    def to_dict(self) -> Dict:
        return {
            "cell_type": self.cell_type,
            "cell_ontology_id": self.cell_ontology_id,
            "reasoning": self.reasoning,
            "expected_markers": self.expected_markers,
            "confidence_estimate": self.confidence_estimate,
            "distinguishing_features": self.distinguishing_features
        }


@dataclass
class EvidenceItem:
    """Single piece of evidence for/against a hypothesis"""
    source: str  # e.g., "cellmarker", "panglaodb", "expression_match"
    content: str
    confidence: float  # 0.0-1.0
    supports: bool  # True=supporting, False=conflicting
    
    def to_dict(self) -> Dict:
        return {
            "source": self.source,
            "content": self.content,
            "confidence": self.confidence,
            "supports": self.supports
        }


@dataclass
class MCTSNode:
    """MCTS tree node for annotation search"""
    # Identity
    node_id: str
    cluster_id: str
    hypothesis: AnnotationHypothesis
    
    # Tree structure
    parent: Optional['MCTSNode'] = None
    children: List['MCTSNode'] = field(default_factory=list)
    
    # MCTS statistics
    visits: int = 0
    total_reward: float = 0.0
    
    # Annotation-specific data
    evidence_chain: List[EvidenceItem] = field(default_factory=list)
    evaluation_result: Optional[Dict] = None
    
    # Cached values
    _ucb_score: Optional[float] = None
    
    @property
    def avg_reward(self) -> float:
        return self.total_reward / self.visits if self.visits > 0 else 0.0
    
    @property
    def confidence(self) -> float:
        """Derived confidence from evaluation result"""
        if self.evaluation_result and "weighted_score" in self.evaluation_result:
            return self.evaluation_result["weighted_score"]
        return self.hypothesis.confidence_estimate
    
    def update(self, reward: float, evidence: List[EvidenceItem], 
               evaluation: Optional[Dict] = None):
        """Update node with search results"""
        self.visits += 1
        self.total_reward += reward
        self.evidence_chain.extend(evidence)
        if evaluation:
            self.evaluation_result = evaluation
        self._ucb_score = None  # Invalidate cache
    
    def ucb_score(self, exploration_weight: float = 1.414, 
                  parent_visits: Optional[int] = None) -> float:
        """Calculate UCB1 score for node selection"""
        if self._ucb_score is not None:
            return self._ucb_score
        
        if parent_visits is None and self.parent:
            parent_visits = self.parent.visits
        elif parent_visits is None:
            parent_visits = 1
        
        if self.visits == 0:
            return float('inf')  # Unvisited nodes get priority
        
        exploitation = self.avg_reward
        exploration = exploration_weight * math.sqrt(
            math.log(parent_visits) / self.visits
        )
        
        self._ucb_score = exploitation + exploration
        return self._ucb_score
    
    def add_child(self, child: 'MCTSNode'):
        child.parent = self
        self.children.append(child)
    
    def is_fully_expanded(self, max_branches: int) -> bool:
        return len(self.children) >= max_branches
    
    def to_summary(self) -> Dict:
        """Lightweight summary for LangGraph state (avoid heavy objects)"""
        return {
            "node_id": self.node_id,
            "cluster_id": self.cluster_id,
            "cell_type": self.hypothesis.cell_type,
            "confidence": self.confidence,
            "visits": self.visits,
            "avg_reward": self.avg_reward,
            "evidence_count": len(self.evidence_chain),
            "child_count": len(self.children)
        }


def select_node(root: MCTSNode, exploration_weight: float = 1.414) -> MCTSNode:
    """Select node for expansion using UCB1"""
    current = root
    
    while current.children:
        # Calculate UCB for all children
        ucb_scores = [
            (child, child.ucb_score(exploration_weight, current.visits))
            for child in current.children
        ]
        # Select child with highest UCB
        current = max(ucb_scores, key=lambda x: x[1])[0]
    
    return current


def backpropagate(node: MCTSNode, reward: float):
    """Propagate reward up the tree"""
    current = node
    while current:
        current.update(reward, [], None)  # Only update stats, no new evidence
        current = current.parent


def extract_best(root: MCTSNode) -> MCTSNode:
    """Extract the node with highest average reward"""
    def traverse(node: MCTSNode) -> MCTSNode:
        best = node
        for child in node.children:
            candidate = traverse(child)
            if candidate.avg_reward > best.avg_reward:
                best = candidate
        return best
    return traverse(root)


def compute_annotation_value(evaluation: Dict, weights: Dict[str, float]) -> float:
    """Compute weighted value from evaluation dimensions"""
    dims = evaluation.get("dimension_scores", {})
    return sum(
        dims.get(dim, 0.0) * weight 
        for dim, weight in weights.items()
    )


def generate_node_id(cluster_id: str, path: List[str]) -> str:
    """Generate unique node ID from cluster and search path"""
    path_hash = hash(tuple(path)) & 0xFFFF
    return f"{cluster_id}_n{path_hash:04x}"

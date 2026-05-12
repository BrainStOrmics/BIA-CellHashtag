# LATS Search Skills

Monte Carlo Tree Search utilities for LATS-based cell annotation in CellHashtag.

## Module Structure

```
skills/lats_search/
├── mcts_core.py          # Core MCTS data structures and algorithms
├── test_mcts_core.py     # Unit tests for mcts_core
└── README.md             # This file
```

## Usage

### Basic MCTS Search

```python
from skills.lats_search.mcts_core import (
    MCTSNode, AnnotationHypothesis, 
    select_node, backpropagate, extract_best
)

# Create root
root = MCTSNode(
    node_id="root",
    cluster_id="cluster_001", 
    hypothesis=AnnotationHypothesis("UNASSIGNED", None, "root", [])
)

# Search loop
for _ in range(n_iterations):
    # Selection
    node = select_node(root, exploration_weight=1.414)
    
    # Expansion (LLM generates hypotheses)
    new_hyps = llm_generate_hypotheses(node.hypothesis)
    for hyp in new_hyps:
        child = MCTSNode(generate_id(), cluster_id, hyp, parent=node)
        node.add_child(child)
    
    # Simulation (validate hypotheses)
    reward = validate_and_score(child.hypothesis)
    
    # Backpropagation
    backpropagate(child, reward)

# Extract best result
best = extract_best(root)
```

### Integration with LangGraph

The `graphs/low_hierarchy/tree_search/` module provides a ready-to-use LangGraph:

```python
from graphs.low_hierarchy.tree_search import create_graph, LATSState

# Build graph
graph = create_graph(config)

# Prepare state
state = LATSState(
    cluster_id="cluster_001",
    cluster_markers=["CD3D", "CD4", "CCR7"],
    tissue_source="lung",
    config_path="config/config.yaml",
    prompts_dir="prompts/low_hierarchy",
    # ... other required fields
)

# Execute
result = graph.invoke(state)
annotation = result["best_annotation"]
```

## Configuration

See `config/config.yaml` for LATS-specific parameters:

```yaml
lats:
  search_params:
    n_iterations: 10
    exploration_weight: 1.414
    max_branches: 3
    confidence_threshold: 0.7
    early_stop_threshold: 0.95
  
  value_function:
    weights:
      marker_match: 0.4
      ontology_consistency: 0.3
      evidence_diversity: 0.2
      llm_confidence: 0.1
```

## Testing

```bash
# Run unit tests
python -m pytest skills/lats_search/test_mcts_core.py -v

# Run with coverage
python -m pytest skills/lats_search/ --cov=skills.lats_search
```

## Design Notes

1. **State Lightweight**: MCTSNode.to_summary() provides a lightweight dict for LangGraph state. Full tree serialization should use external storage for large searches.

2. **Parallel Validation**: The `validate_hypotheses` node is designed for parallel execution. In production, use LangGraph's `Send` API to fan out validation tasks.

3. **Caching**: Evidence results are cached in `evidence_cache` to avoid redundant API calls to CellWiki or other databases.

4. **Early Termination**: The search can terminate early if confidence exceeds `early_stop_threshold`, saving compute resources.

## References

- Zhou, A., et al. (2024). Language Agent Tree Search Unifies Reasoning, Acting, and Planning. ICML.
- [LangGraph Examples: LATS](https://github.com/langchain-ai/langgraph/tree/main/examples/lats)
- [CellHashtag Architecture](../../ARCHITECTURE.md)

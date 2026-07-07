---
name: lats-search
description: LATS (Look-Ahead Tree Search) for cell type annotation using MCTS. Explores hypothesis space via iterative expansion, validation, and backpropagation.
allowed-tools: Read, Bash
compatibility: Single-cell cluster annotation with uncertain cell types
---

# LATS MCTS Protocol

## When to Use

For each cluster that needs annotation. Run MCTS to find the best cell type hypothesis by iteratively expanding, validating, and scoring candidates.

## Algorithm

### 1. Initialize
Create root `MCTSNode` with hypothesis `{"cell_type": "UNASSIGNED", "expected_markers": cluster_markers}`.

### 2. Select (UCB1)
Traverse from root to leaf, at each node pick the child with highest UCB1 score:
```
UCB1 = avg_reward + exploration_weight * sqrt(log(parent_visits) / visits)
```
exploration_weight = 1.414 (sqrt(2)). Unvisited nodes get `inf` priority.

### 3. Expand (LLM)
Generate 1-3 new hypotheses using the LLM:
- Consider tissue context and experimental conditions
- Ensure biological plausibility
- Diversify candidates from existing hypotheses
- Each hypothesis: `{"cell_type": str, "expected_markers": list[str], "reasoning": str, "confidence_estimate": float}`

### 4. Validate (CellWiki)
For each hypothesis, check expected markers against CellWiki:
- Query known cell types for overlapping markers
- Compute Jaccard similarity between observed and expected markers
- Cache evidence in state

### 5. Evaluate (LLM)
Score each hypothesis on 4 dimensions:
- marker_match (0.4 weight): overlap between cluster markers and expected markers
- ontology_consistency (0.3): does cell type fit tissue context?
- evidence_diversity (0.2): multiple independent evidence sources
- llm_confidence (0.1): LLM self-assessment

### 6. Backpropagate
Update visit counts and total reward from evaluated node back to root.

### 7. Terminate
Stop when:
- Best confidence >= early_stop_threshold (0.95): confident
- Iterations >= max_iterations (10): exhausted
- Best confidence >= confidence_threshold (0.7): acceptable

### 8. Extract
Return the node with highest `avg_reward` as best annotation:
```json
{"cell_type": str, "confidence": float, "reasoning": str, "evidence": [...]}
```

## Configuration

Read from `config.yaml` under `lats.search_params` and `lats.value_weights`.

## Files

- `core/mcts.py`: `MCTSNode`, `select_node`, `backpropagate`, `extract_best`
- `skills/cellwiki/`: `query_markers_for_genes`, `query_cell_type`

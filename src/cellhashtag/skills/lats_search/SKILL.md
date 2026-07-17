---
name: lats-search
description: DAG-MCTS for cell type annotation. Search space is a DAG (Merge + ReEmbed create multi-parent nodes). Uses P-UCT selection, Soft-Bellman backprop, confidence discounting, and transposition table. Per-iteration tree viz emitted as interactive HTML.
allowed-tools: Read, Bash
compatibility: Single-cell cluster annotation with uncertain cell types
---

# DAG-MCTS Protocol (v4.0+)

Replaces the tree-only UCB1 MCTS in `core/mcts.py` (deleted). Search space is a DAG: cross-branch `Merge` and `ReEmbed` create multi-parent nodes. Transposition table dedups convergent states.

## When to Use

For each cluster needing annotation. Run DAG-MCTS to find the best cell type hypothesis via:
1. P-UCT selection with FPU + priority heuristic (Relevance × Heterogeneity × Uncertainty)
2. LLM-driven expansion into the DAG action space
3. Adversarial LLM evaluation with 5-component reward + sliding-window normalization
4. Soft-Bellman backprop with temperature τ-scheduling
5. Confidence discounting via DAG in-degree: λ_shared = 1/(1+α·(|Parents|-1))
6. Reflection schema (diagnosis / assessment / suggested_next / hallucination_risk)
7. Per-iteration PyVis HTML tree viz dump

## Cold Start (3-stage)

Before the main loop, `cold_start()` builds the root node:
1. **Metadata match** — tissue/organ hints → pick prior template (pbmc/tumor/brain/...)
2. **Quick profiling** — scanpy Leiden to estimate cluster heterogeneity
3. **LLM prior injection** — suggest expected lineages and marker panels

Output: root `DAGNode` with informed `semantic_state` and `topology_state`.

## Algorithm

### 1. Initialize
Create root `DAGNode` via `cold_start()` with `state_hash = sha256(topology | semantic)`.

### 2. Select (P-UCT + FPU)
Traverse from root to leaf. At each node pick child with highest P-UCT score:
```
P-UCT(c) = FPU(c) + cpuct * prior_P(c) * sqrt(N_parent) / (1 + N_c)
```
FPU for unvisited child with priority score `p`:
- `p >= theta_high (0.7)` → `V_root + delta`
- `p < theta_low (0.3)` → `V_root - delta`
- else → `V_root`

### 3. Expand (LLM action shortlist)
LLM proposes 1-3 actions from `{Split, Merge, Freeze, AssignLabel}` with weights. Prior_P = (weight / sum) × (1-ε) + ε/k.

### 4. Transposition
On expand: if `state_hash(child)` exists in `TranspositionTable`, reuse the node and add new parent (creating DAG edge). Otherwise create new node.

### 5. Evaluate (LLM adversarial + reflection)
5-component reward:
- purity (0.25) — cluster internal homogeneity
- specificity (0.20) — marker exclusivity vs other clusters
- context (0.20) — tissue ontology consistency
- lats (0.20) — reflection self-assessment
- known_marker (0.15) — CellWiki canonical overlap

Plus overfragmentation penalty. Scores normalized via `SlidingWindowNormalizer` (window=50). Potential shaping (Ng) applied.

LLM also emits reflection JSON:
```json
{"diagnosis": "...", "assessment": 0.0-1.0, "suggested_next": "...",
 "hallucination_risk": "low|medium|high", "should_continue": true}
```

### 6. Backpropagate (Soft-Bellman)
For each node on path (leaf → root), with depth d:
```
τ_d = τ0 / (1 + γ·d)
V(node) = τ_d · log(Σ_c w_c · exp(Q(c)/τ_d))   [Log-Sum-Exp]
λ_shared = 1 / (1 + α·(|Parents|-1))
λ_halluc = risk_weight(hallucination_risk) ∈ {1.0, 0.6, 0.2}
Δ = r_immediate + λ_shared·λ_halluc·V_successor − Q(node)
Q(node) += Δ / max(visits, 1)
visits += 1
```

### 7. Tree Viz
After each iteration, emit `output/trees/<cluster_id>/iter_{n:03d}.html` via PyVis + NetworkX:
- Node color = hallucination_risk (green/yellow/red)
- Node size = 10 + 8·log1p(visits)
- Edge width = 0.5 + 3·prior_P
- Hierarchical UD layout

### 8. Terminate
Stop when:
- `best_score >= early_stop_threshold (0.95)`: confident
- `iterations >= max_iterations (20)`: exhausted

### 9. Extract
Return the `DAGNode` with highest `q_value` (filtered to `visits > 0`).

## Configuration

Read from `config.yaml` under `dag_mcts:` section. All params exposed:
- `max_iterations`, `cpuct`, `tau0`, `tau_gamma`, `alpha`, `epsilon`
- `early_stop_threshold`, `confidence_threshold`
- `fpu_delta`, `theta_high`, `theta_low`
- `reward_weights` (5-component dict)
- `enable_tree_viz` (bool)

## Files

- `core/dag_mcts.py`: `DAGNode`, `TranspositionTable`, `run_dag_mcts`, P-UCT, Soft-Bellman
- `core/reflection.py`: `ReflectionSchema`, `parse_reflection`
- `core/reward.py`: `RewardBreakdown`, `SlidingWindowNormalizer`, `compute_reward`
- `core/cold_start.py`: `cold_start()`, `ColdStartResult`
- `core/tree_viz.py`: `TreeVisualizer` (PyVis HTML per iteration)
- `core/lats_loop.py`: `run_lats_search()` — wires DAG-MCTS to DeepAgents LLM calls
- `skills/cellwiki/`: `query_markers_for_genes`, `query_cell_type`

## Dual-Track Architecture

- **Topology track**: bio tools (scanpy, Wilcoxon marker test, CellWiki queries)
- **Semantic track**: LATS-style reflection with hallucination risk flagging

The two tracks compose at `evaluate()` — bio scores feed reward, reflection feeds risk discount.

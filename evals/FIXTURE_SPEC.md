# Fixture Specification

Every fixture is JSONL (one JSON record per line). Each record has the same
envelope regardless of category:

```json
{
  "id": "ann_pbmc_cd4_naive_01",
  "category": "annotation",
  "version": "2026-06-11",
  "input": { ... },
  "expected": { ... },
  "metadata": { ... }
}
```

## Envelope Fields

| Field | Type | Required | Notes |
|---|---|---|---|
| `id` | str | yes | Globally unique. `<category>_<dataset>_<case>_<nn>` |
| `category` | enum | yes | `cluster_quality` \| `annotation` \| `mcts` \| `tool_use` \| `ca_aa` |
| `version` | ISO date | yes | Fixture schema version; bump on breaking change |
| `input` | object | yes | Everything the system under test receives |
| `expected` | object | yes | Gold output the scorer compares against |
| `metadata` | object | no | Provenance, difficulty, tissue, species |

## Per-Category Schemas

### `cluster_quality`

Input: a pickled AnnData path + clustering parameters.
Expected: acceptable ranges for ASW, modularity, cluster count.

```json
{
  "id": "cq_pbmc3k_res06_01",
  "category": "cluster_quality",
  "version": "2026-06-11",
  "input": {
    "adata_pickle": "evals/fixtures/cluster_quality/data/pbmc3k.pkl",
    "resolution": 0.6,
    "n_neighbors": 15,
    "n_pcs": 30
  },
  "expected": {
    "asw_min": 0.25,
    "asw_max": 0.70,
    "modularity_min": 0.30,
    "cluster_count_range": [8, 16],
    "min_cells_per_cluster": 20
  },
  "metadata": {"dataset": "pbmc3k", "species": "human", "tissue": "PBMC"}
}
```

Provenance: 10x Genomics pbmc3k, Seurat tutorial canonical clustering.

### `annotation`

Input: cluster top markers + tissue context.
Expected: gold cell type + canonical tier-1 markers.

```json
{
  "id": "ann_pbmc_cd4_naive_01",
  "category": "annotation",
  "version": "2026-06-11",
  "input": {
    "cluster_id": "c_7",
    "top_markers": ["CD3D","CD3E","CD4","CCR7","LEF1","SELL","IL7R","TCF7","MAL","LCK"],
    "tissue": "PBMC",
    "species": "human",
    "mean_expression_top5": [2.1, 1.9, 1.4, 0.8, 0.7]
  },
  "expected": {
    "cell_type": "CD4 Naive T",
    "lineage": "T cell",
    "tier1_markers": ["CD3D","CD3E","CD4","CCR7","LEF1"],
    "acceptable_alternates": ["CD4 T naive", "Naive CD4+ T"],
    "excluded_types": ["CD8 T", "B cell", "NK", "Monocyte"]
  },
  "metadata": {"source": "Zheng2017", "n_cells": 412, "difficulty": "easy"}
}
```

Provenance: Zheng et al. 2017 sorted PBMC; CellMarker 2.0 canonical lists.

### `mcts`

Input: a deterministic MCTS root state (serialized MCTSNode tree).
Expected: best action at root under known reward landscape.

```json
{
  "id": "mcts_cd4_naive_converge_01",
  "category": "mcts",
  "version": "2026-06-11",
  "input": {
    "root": {"node_id":"r","hypothesis":"UNASSIGNED","visits":0,"reward":0.0,"children":["a","b","c"]},
    "nodes": [
      {"node_id":"a","hypothesis":"CD4 Naive T","visits":5,"reward":3.8,"children":["a1"]},
      {"node_id":"b","hypothesis":"CD4 Memory T","visits":3,"reward":1.2,"children":[]},
      {"node_id":"c","hypothesis":"CD8 T","visits":1,"reward":0.1,"children":[]},
      {"node_id":"a1","hypothesis":"CD4 Naive T (sub)","visits":2,"reward":1.5,"children":[]}
    ],
    "ucb_c": 1.41,
    "gamma": 0.9
  },
  "expected": {
    "best_action": "a",
    "selection_path": ["r","a"],
    "best_value_min": 0.70,
    "max_depth_effective": 3,
    "exploration_share_min": 0.20
  },
  "metadata": {"scenario": "clear_winner", "difficulty": "easy"}
}
```

Provenance: synthetic, hand-computed UCB1 values.

### `tool_use`

Input: a natural-language or structured request the agent must satisfy by calling a tool.
Expected: which tool, what args.

```json
{
  "id": "tu_cellwiki_cd4_01",
  "category": "tool_use",
  "version": "2026-06-11",
  "input": {
    "task": "validate_markers",
    "cell_type_hypothesis": "CD4 Naive T",
    "observed_markers": ["CD3D","CD3E","CD4","CCR7","LEF1"],
    "tissue": "PBMC"
  },
  "expected": {
    "tool": "cellwiki.query_markers_for_genes",
    "required_args": {"genes": ["CD3D","CD3E","CD4","CCR7","LEF1"]},
    "optional_args": {"tissue": "PBMC"},
    "must_not_call": ["cellwiki.list_cell_types"]
  },
  "metadata": {"difficulty": "easy"}
}
```

```json
{
  "id": "tu_clustering_quality_sil_01",
  "category": "tool_use",
  "version": "2026-06-11",
  "input": {
    "task": "assess_clustering",
    "adata_pickle": "evals/fixtures/tool_use/data/small.pkl"
  },
  "expected": {
    "tools_called_set": [
      "clustering_quality.silhouette_check",
      "clustering_quality.modularity_check"
    ],
    "min_tools_called": 2
  },
  "metadata": {}
}
```

### `ca_aa`

Input: a snapshot of the CA-AA loop state.
Expected: the next action the orchestrator should take.

```json
{
  "id": "caaa_merge_same_type_01",
  "category": "ca_aa",
  "version": "2026-06-11",
  "input": {
    "alternation_count": 1,
    "max_alternations": 3,
    "clusters": [
      {"id":"c_0","cell_type":"CD4 Naive T","confidence":0.91,"n_cells":120},
      {"id":"c_3","cell_type":"CD4 Naive T","confidence":0.88,"n_cells":95},
      {"id":"c_5","cell_type":"B cell","confidence":0.82,"n_cells":210},
      {"id":"c_9","cell_type":"Unknown","confidence":0.22,"n_cells":40}
    ]
  },
  "expected": {
    "action": "merge",
    "merge_candidates": [["c_0","c_3"]],
    "subcluster_targets": ["c_9"],
    "converged": false
  },
  "metadata": {"scenario": "merge+subcluster"}
}
```

```json
{
  "id": "caaa_converged_01",
  "category": "ca_aa",
  "version": "2026-06-11",
  "input": {
    "alternation_count": 2,
    "max_alternations": 3,
    "clusters": [
      {"id":"c_0","cell_type":"CD4 Naive T","confidence":0.85},
      {"id":"c_3","cell_type":"B cell","confidence":0.90},
      {"id":"c_5","cell_type":"NK","confidence":0.78}
    ]
  },
  "expected": {
    "action": "converge",
    "merge_candidates": [],
    "subcluster_targets": [],
    "converged": true
  },
  "metadata": {"scenario": "all_high_conf"}
}
```

## Fixture Provenance Rules

- Every fixture must declare its gold source in `metadata.source`
- Synthetic fixtures are marked `source: "synthetic"` and require a hand-computed
  expected output committed alongside the derivation in `fixtures/<category>/derivations/`
- Human-labeled fixtures require ≥2 annotators; disagreements resolved by a senior
  biologist and recorded in `fixtures/<category>/labels.csv`
- Minimum fixture counts per category for CI:
  - `cluster_quality`: 5
  - `annotation`: 50 (10 easy / 30 medium / 10 hard)
  - `mcts`: 20
  - `tool_use`: 15
  - `ca_aa`: 15

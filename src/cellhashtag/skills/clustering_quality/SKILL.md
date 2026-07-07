---
name: clustering-quality
description: Assess single-cell clustering quality using silhouette score, modularity, and phiclust. Returns recommendations to adjust resolution or accept current clustering.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Clustering Quality Assessment

## When to Use

After running Leiden or Louvain clustering on single-cell data. Use before cell type annotation to ensure cluster boundaries are biologically meaningful.

## Tools

### `silhouette_check(adata, cluster_key, use_rep)`

Computes Average Silhouette Width (ASW) on the PCA embedding.

- ASW > 0.5: clear cluster structure
- ASW 0.25-0.5: moderate
- ASW < 0.25: recommend adjusting resolution

Returns: `{"tool": "silhouette", "asw": float, "per_cluster": dict, "recommendation": str, "details": str}`

### `modularity_check(adata, cluster_key)`

Computes graph modularity on the connectivities graph.

- Modularity > 0.3: meaningful community structure
- Modularity > 0.5: strong community structure
- Modularity < 0.2: recommend adjusting resolution

Requires `adata.obsp["connectivities"]` from `sc.pp.neighbors()`.

Returns: `{"tool": "modularity", "modularity": float, "recommendation": str, "details": str}`

### `phiclust_check(adata, cluster_key)`

R-based phi score analysis for hidden substructure. Requires `rpy2`.

- Phi > 0 for >50% of clusters: increase resolution
- No clusters with phi > 0: possible over-clustering

## Decision Logic

Majority vote across available tools:
- If majority says "adjust": multiply resolution by 1.5 (if under-clustered) or 0.7 (if over-clustered)
- If majority says "ok": accept current clustering
- Resolution clamped to [0.05, 2.0]

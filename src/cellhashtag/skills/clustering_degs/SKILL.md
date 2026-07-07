---
name: clustering-degs
description: Find differentially expressed genes (markers) for clusters using wilcoxon, t-test, or logistic regression. Includes consensus markers and resolution scanning.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Clustering Differential Expression

## When to Use

After clustering to identify marker genes for each cluster, or before annotation to refine cluster resolution.

## Functions

### `find_markers(adata, groupby, method="wilcoxon", top_n=20) -> dict`

Finds DEGs per cluster.

Parameters:
- `method`: "wilcoxon" (default), "t-test", or "logreg"
- `top_n`: number of top markers per cluster

Returns dict keyed by cluster ID with list of gene names.

### `consensus_markers(adata, groupby, methods=("wilcoxon", "t-test"), min_methods=2, top_n=20) -> dict`

Finds markers that appear in results from multiple statistical methods.

Parameters:
- `methods`: tuple of methods to compare
- `min_methods`: minimum number of methods that must agree (default 2)

Returns dict of consensus markers per cluster.

### `resolution_scan(adata, resolutions, cluster_key="leiden", quality_fn="silhouette") -> dict`

Scans multiple resolution values and returns quality scores.

Returns `{resolution: score}` dict and recommends best resolution.

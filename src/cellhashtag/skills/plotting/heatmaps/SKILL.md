---
name: plotting-heatmaps
description: Generate marker gene heatmaps with hierarchical clustering dendrograms.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Heatmap Plots

## When to Use

To visualize marker gene expression patterns across clusters for annotation review.

## Functions

### `heatmap(adata, genes, groupby, output_dir=None, dendrogram=True, show_values=False) -> str`

Creates marker gene heatmap with optional dendrogram.
Saves to `{output_dir}/heatmap_{groupby}.png`

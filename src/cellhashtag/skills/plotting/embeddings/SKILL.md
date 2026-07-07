---
name: plotting-embeddings
description: Generate publication-quality UMAP/t-SNE embedding plots with Okabe-Ito color palette. Outputs PDF+PNG.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq, spatial
---

# Embedding Plots

## When to Use

After clustering or annotation to visualize cell distributions.

## Functions

### `plot_embedding(adata, embedding_key, color_by=None, output_dir=None, title=None, colors=None) -> str`

Generates publication-quality embedding plot.

Parameters:
- `embedding_key`: "X_umap", "X_tsne", or any key in adata.obsm
- `color_by`: Column name in adata.obs for coloring (e.g. "leiden", "Cell#")
- `output_dir`: Directory to save files. If None, returns without saving
- `title`: Plot title. Auto-generated if None
- `colors`: Optional list of specific colors (uses Okabe-Ito palette by default)

Saves as both PDF and PNG in output_dir.
Uses Okabe-Ito colorblind-friendly palette.
DPI=300 for PNG.

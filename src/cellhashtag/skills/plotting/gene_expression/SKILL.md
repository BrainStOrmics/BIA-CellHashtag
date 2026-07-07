---
name: plotting-gene-expression
description: Generate violin plots and dot plots for gene expression visualization across cell clusters.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Gene Expression Plots

## When to Use

During marker gene validation and cell type annotation review.

## Functions

### `violin_plot(adata, genes, groupby, output_dir=None) -> str`

Creates violin plot showing expression distribution of specified genes across groups.
Saves to `{output_dir}/violin_{groupby}.png`

### `dot_plot(adata, genes, groupby, output_dir=None) -> str`

Creates dot plot showing average expression and fraction of cells expressing each gene per group.
Dot size = fraction of cells expressing gene.
Dot color = average expression level.
Saves to `{output_dir}/dotplot_{groupby}.png`

---
name: quality-control
description: Compute QC metrics, filter cells, detect doublets for single-cell RNA-seq data.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Quality Control

## When to Use

Before clustering to remove low-quality cells and doublets.

## Functions

### `compute_qc_metrics(adata) -> dict`

Computes standard QC metrics:
- `n_genes`: number of detected genes per cell
- `total_counts`: total UMI counts per cell
- `pct_mito`: percentage of mitochondrial genes
- `pct_ribo`: percentage of ribosomal genes

Returns metrics added to `adata.obs`: `n_genes`, `total_counts`, `pct_mito`, `pct_ribo`

Mitochondrial genes detected by prefix "MT-" (human) or "mt-" (mouse).
Ribo genes by prefixes "RPS", "RPL" (ribo protein) and "RPSA".

### `filter_cells(adata, min_genes=200, max_genes=5000, min_counts=500, max_counts=50000, max_mito_pct=20) -> AnnData`

Filters cells based on QC thresholds. Returns subsetted AnnData.
- Cells with too few genes = likely empty droplets
- Cells with too many genes = likely doublets
- High mito percentage = dying cells

### `detect_doublets(adata, expected_cells=None) -> list[int]`

Heuristic doublet detection based on gene count outliers.
Cells with `n_genes > Q3 + 1.5*IQR` flagged as potential doublets.
Returns list of cell indices flagged as doublets.

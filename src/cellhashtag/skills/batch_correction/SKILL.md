---
name: batch-correction
description: Correct batch effects in single-cell data using Harmony, BBKNN, or ComBat.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Batch Correction

## When to Use

When integrating datasets from multiple batches, donors, or experimental runs.

## Functions

### `correct_harmony(adata, batch_key, use_rep="X_pca") -> AnnData`

Runs Harmony batch correction on PCA embedding.
Returns AnnData with corrected embedding in `adata.obsm["X_pca_harmony"]`.
Requires `harmonypy` package.

### `correct_bbknn(adata, batch_key, neighbors_within_batch=15) -> AnnData`

Runs BBKNN to build batch-aware neighbor graph.
Modifies `adata.obsp["connectivities"]` and `adata.obsp["distances"]`.
Requires `bbknn` package.

### `correct_combat(adata, batch_key, key_added="X_pca_combat") -> AnnData`

Runs ComBat batch correction.
Returns AnnData with corrected data in `adata.obsm[key_added]`.
Uses `scanpy.external.pp.combat`.

---
name: spatial-loader
description: Load spatial transcriptomics data from 10x Visium, Nanostring MERFISH, BGI Stereo-seq, and generic spatial formats.
allowed-tools: Read, Bash
compatibility: spatial transcriptomics
---

# Spatial Data Loaders

## When to Use

To load spatial transcriptomics data into AnnData format for downstream analysis.

## Functions

### `load_visium(path, sample_id=None) -> AnnData`

Loads 10x Visium data from the output directory.
Expects `filtered_feature_bc_matrix/` and `spatial/` subdirectories.
Uses `scanpy.read_visium()`.
Adds `spatial` to `adata.obsm`.

### `load_merfish(path) -> AnnData`

Loads Nanostring MERFISH data from CSV/parquet.
Expects gene expression matrix with spatial coordinates (x, y).
Returns AnnData with coordinates in `adata.obsm["spatial"]`.

### `load_stereoseq(path) -> AnnData`

Loads BGI Stereo-seq data from GEM file or sparse matrix.
Returns AnnData with coordinates in `adata.obsm["spatial"]`.

### `load_spatial(expression_path, coordinates_path=None, spatial_coords=None) -> AnnData`

Generic spatial data loader.
Accepts expression matrix (CSV/parquet) and either a coordinates file or numpy array of (x, y) positions.

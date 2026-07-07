---
name: format-conversion
description: Load single-cell RNA-seq data from h5ad, rds, or 10x mtx formats into AnnData objects.
allowed-tools: Read, Bash
compatibility: h5ad, rds, sce, 10x mtx
---

# Data Format Conversion

## Supported Formats

### h5ad (`load_h5ad(path)`)
Native AnnData format. Direct load via `scanpy.read_h5ad()`.
- Preserves `.X`, `.obs`, `.var`, `.obsm`, `.uns`, `.obsp`
- Fastest format, recommended

### 10x mtx (`load_mtx(path)`)
10x Genomics matrix directory. Expects:
- `matrix.mtx` (or `matrix.mtx.gz`)
- `genes.tsv` (or `features.tsv`)
- `barcodes.tsv`
Loaded via `scanpy.read_10x_mtx(path, var_names="gene_symbols")`.

### rds (`load_rds(path)`)
R Seurat object saved as .rds. Requires `rpy2` and `anndata2ri` or `zellkonverter`.
- Converts Seurat assay to AnnData format
- Most complex conversion, fallback may be needed

### sce (`load_sce(path)`)
SingleCellExperiment object saved as .rds. Requires `zellkonverter` and `rpy2`.
Loaded via `zellkonverter.readRDS(path)`.
- Converts SCE assay, rowData, colData to AnnData equivalents
- More reliable than Seurat→AnnData conversion for SCE-native data

## Return Value
All loaders return an `AnnData` object ready for downstream analysis.

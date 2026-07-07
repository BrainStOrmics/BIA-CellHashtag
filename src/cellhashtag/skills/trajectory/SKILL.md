---
name: trajectory-inference
description: Infer cell trajectories using PAGA graph-based trajectory inference.
allowed-tools: Read, Bash
compatibility: scRNA-seq, snRNA-seq
---

# Trajectory Inference

## When to Use

To understand developmental trajectories and cell state transitions within a dataset.

## Functions

### `paga_trajectory(adata, cluster_key="leiden", plot=True, output_dir=None) -> dict`

Runs PAGA trajectory inference on AnnData.

Requires `sc.pp.neighbors()` and `sc.tl.leiden()` to have been run.

Returns dict with:
- `connectivities_tree`: PAGA tree as networkx graph
- `plot_path`: path to saved PAGA plot (if output_dir provided)

Parameters:
- `cluster_key`: Column in adata.obs for clustering
- `plot`: Whether to generate PAGA plot
- `output_dir`: Directory to save plot

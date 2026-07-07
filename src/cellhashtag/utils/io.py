"""AnnData I/O, expression summary, and marker extraction utilities."""

import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def load_adata(path: str, fmt: str = "h5ad") -> AnnData:
    loaders = {"h5ad": sc.read_h5ad, "mtx": _load_mtx}
    if fmt not in loaders:
        raise ValueError(f"Unsupported format: {fmt}. Use 'h5ad' or 'mtx'.")
    return loaders[fmt](path)


def _load_mtx(path: str) -> AnnData:
    return sc.read_10x_mtx(path, var_names="gene_symbols")


def save_adata(adata: AnnData, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)


def adata_to_pickle(adata: AnnData, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(adata, f)


def adata_from_pickle(path: str) -> AnnData:
    with open(path, "rb") as f:
        return pickle.load(f)


def extract_markers(adata: AnnData, cluster_key: str, cluster_id: str, top_n: int = 20) -> list[str]:
    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, reference="rest")
    raw = list(adata.uns["rank_genes_groups"]["names"][cluster_id][:top_n])
    return [g for g in raw if g and not str(g).startswith("MT-")]


def expression_summary(adata: AnnData, genes: list[str], exp_cutoff: float = 0.01) -> str:
    lines = ["==========\n", "| gene name | expression level | expression ratio (%) |\n", "|-----------|------------------|----------------------|\n"]
    for gene in genes:
        if gene not in adata.var_names:
            continue
        mat = adata[:, adata.var_names == gene].X
        if mat.shape[1] != 1:
            continue
        if hasattr(mat, "toarray"):
            mat = mat.toarray()
        level = float(np.mean(mat))
        ratio = adata[mat > exp_cutoff].n_obs / adata.n_obs * 100
        lines.append(f"| {gene} | {level:.2f} | {ratio:.2f} |\n")
    lines.append("==========\n")
    return "".join(lines)


def infer_omics_type(adata: AnnData) -> str:
    obs_cols = set(adata.obs.columns)
    obsm_keys = set(adata.obsm.keys())
    if any("spatial" in k for k in obsm_keys) or "imagerow" in obs_cols:
        return "spRNA"
    if any(k in obs_cols for k in ["protein_counts", "adt_counts"]):
        return "CITE"
    if any("peak" in k for k in obsm_keys):
        return "ATAC"
    return "scRNA"


def perceive_data(adata: AnnData) -> dict:
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "gene_names": list(adata.var_names[:10]),
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "is_sparse": hasattr(adata.X, "toarray"),
        "omics_type": infer_omics_type(adata),
    }


def df2markdown_table(df: pd.DataFrame) -> str:
    header = "| " + " | ".join(df.columns) + " |"
    separator = "| " + " | ".join(["---"] * len(df.columns)) + " |"
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for _, row in df.iterrows()]
    return header + "\n" + separator + "\n" + "\n".join(rows)

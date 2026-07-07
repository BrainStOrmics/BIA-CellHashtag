"""Nanostring MERFISH spatial data loader."""
from typing import Any
import os
import pandas as pd
import numpy as np
import anndata as ad


def load_merfish(path: str) -> Any:
    if not os.path.exists(path):
        raise FileNotFoundError(f"MERFISH file not found: {path}")

    if path.endswith(".parquet"):
        df = pd.read_parquet(path)
    elif path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        raise ValueError("MERFISH file must be .csv or .parquet")

    gene_cols = [c for c in df.columns if c not in ("x", "y", "X_center", "Y_center", "cell_id", "barcode")]
    expr = df[gene_cols].values.astype(np.float32)

    coord_keys = [k for k in ("x", "y", "X_center", "Y_center") if k in df.columns]
    coords = df[coord_keys[:2]].values.astype(np.float32) if len(coord_keys) >= 2 else np.zeros((len(df), 2))

    obs = pd.DataFrame(index=df.get("cell_id", [f"cell_{i}" for i in range(len(df))]))
    for k in ("x", "y", "X_center", "Y_center", "barcode"):
        if k in df.columns:
            obs[k] = df[k].values

    adata = ad.AnnData(X=expr, obs=obs)
    adata.var_names = gene_cols
    adata.obsm["spatial"] = coords

    return adata

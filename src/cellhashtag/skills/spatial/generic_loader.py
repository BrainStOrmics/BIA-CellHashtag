"""Generic spatial data loader."""
from typing import Any
import os
import pandas as pd
import numpy as np
import anndata as ad


def load_spatial(
    expression_path: str,
    coordinates_path: str | None = None,
    spatial_coords: np.ndarray | None = None,
) -> Any:
    if not os.path.exists(expression_path):
        raise FileNotFoundError(f"Expression file not found: {expression_path}")

    if expression_path.endswith(".csv"):
        expr_df = pd.read_csv(expression_path, index_col=0)
    elif expression_path.endswith(".parquet"):
        expr_df = pd.read_parquet(expression_path)
    else:
        raise ValueError("Expression file must be .csv or .parquet")

    expr = expr_df.values.astype(np.float32)
    adata = ad.AnnData(X=expr)
    adata.obs_names = expr_df.index.astype(str)
    adata.var_names = expr_df.columns.astype(str)

    if spatial_coords is not None:
        coords = np.asarray(spatial_coords, dtype=np.float32)
    elif coordinates_path:
        if coordinates_path.endswith(".csv"):
            coords_df = pd.read_csv(coordinates_path)
        else:
            coords_df = pd.read_csv(coordinates_path, sep="\t")
        coords = coords_df[["x", "y"]].values.astype(np.float32)
    else:
        raise ValueError("Provide either coordinates_path or spatial_coords")

    adata.obsm["spatial"] = coords
    return adata

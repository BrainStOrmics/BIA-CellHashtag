"""BGI Stereo-seq spatial data loader."""
from typing import Any
import os
import pandas as pd
import numpy as np
import anndata as ad


def load_stereoseq(path: str) -> Any:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Stereo-seq file not found: {path}")

    if path.endswith(".csv"):
        df = pd.read_csv(path)
    elif path.endswith(".tsv"):
        df = pd.read_csv(path, sep="\t")
    elif path.endswith(".parquet"):
        df = pd.read_parquet(path)
    else:
        raise ValueError("Stereo-seq file must be .csv, .tsv, or .parquet")

    required = {"gene", "x", "y"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    adata = ad.AnnData(
        X=df.groupby(["gene", "x", "y"]).size().unstack(fill_value=0).T.values,
    )
    adata.obs_names = df["x"].unique()
    adata.var_names = df["gene"].unique()
    adata.obsm["spatial"] = np.column_stack([df["x"].unique(), df["y"].unique()])[:len(adata)]

    return adata

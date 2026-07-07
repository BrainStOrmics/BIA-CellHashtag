"""10x Visium spatial data loader."""
from typing import Any


def load_visium(path: str, sample_id: str | None = None) -> Any:
    import os
    import scanpy as sc

    if not os.path.exists(path):
        raise FileNotFoundError(f"Visium directory not found: {path}")

    adata = sc.read_visium(path)
    adata.var_names_make_unique()

    return adata

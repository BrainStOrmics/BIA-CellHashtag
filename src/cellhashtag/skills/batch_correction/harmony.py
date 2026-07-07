"""Harmony batch correction."""
from typing import Any


def correct_harmony(adata: Any, batch_key: str, use_rep: str = "X_pca") -> Any:
    try:
        import harmonypy as hm
    except ImportError as exc:
        raise ImportError("harmonypy required: pip install harmonypy") from exc

    import numpy as np

    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    if use_rep not in adata.obsm:
        raise ValueError(f"{use_rep} not found in adata.obsm")

    meta = adata.obs[[batch_key]].copy()
    vars_use = [batch_key]

    ho = hm.run_harmony(adata.obsm[use_rep], meta, vars_use)
    adata.obsm["X_pca_harmony"] = np.array(ho.Z_corr).T

    return adata

"""BBKNN batch correction."""
from typing import Any


def correct_bbknn(adata: Any, batch_key: str, neighbors_within_batch: int = 15) -> Any:
    try:
        import bbknn
    except ImportError as exc:
        raise ImportError("bbknn required: pip install bbknn") from exc

    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    bbknn.bbknn(adata, batch_key=batch_key, neighbors_within_batch=neighbors_within_batch)
    return adata

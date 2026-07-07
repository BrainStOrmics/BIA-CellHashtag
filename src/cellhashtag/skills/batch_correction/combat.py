"""ComBat batch correction."""
from typing import Any


def correct_combat(adata: Any, batch_key: str, key_added: str = "X_pca_combat") -> Any:
    import scanpy as sc
    import numpy as np

    if batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key '{batch_key}' not in adata.obs")

    adata_tmp = adata.copy()
    adata_tmp.X = adata.X.copy() if hasattr(adata.X, "copy") else np.array(adata.X)

    sc.external.pp.combat(adata_tmp, key=batch_key)
    adata.obsm[key_added] = adata_tmp.X

    return adata

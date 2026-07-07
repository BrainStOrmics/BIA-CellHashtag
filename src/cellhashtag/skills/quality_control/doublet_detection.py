import numpy as np

from anndata import AnnData

from .qc_metrics import compute_qc_metrics


def detect_doublets(adata: AnnData, expected_cells: int | None = None) -> list[int]:
    if "n_genes" not in adata.obs.columns:
        compute_qc_metrics(adata)

    n_genes = adata.obs["n_genes"].values

    q1 = np.percentile(n_genes, 25)
    q3 = np.percentile(n_genes, 75)
    iqr = q3 - q1

    threshold = q3 + 1.5 * iqr

    doublet_mask = n_genes > threshold
    doublet_indices = np.where(doublet_mask)[0].tolist()

    return doublet_indices

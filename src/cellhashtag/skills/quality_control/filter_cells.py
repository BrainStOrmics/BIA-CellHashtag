import logging

from anndata import AnnData

from .qc_metrics import compute_qc_metrics

logger = logging.getLogger(__name__)


def filter_cells(
    adata: AnnData,
    min_genes: int = 200,
    max_genes: int = 5000,
    min_counts: int = 500,
    max_counts: int = 50000,
    max_mito_pct: float = 20
) -> AnnData:
    if "n_genes" not in adata.obs.columns:
        compute_qc_metrics(adata)

    initial_count = adata.n_obs

    gene_mask = (adata.obs["n_genes"] >= min_genes) & (adata.obs["n_genes"] <= max_genes)
    count_mask = (adata.obs["total_counts"] >= min_counts) & (adata.obs["total_counts"] <= max_counts)
    mito_mask = adata.obs["pct_mito"] <= max_mito_pct

    combined_mask = gene_mask & count_mask & mito_mask
    filtered_adata = adata[combined_mask].copy()

    removed = initial_count - filtered_adata.n_obs
    logger.info(
        f"Filtered {removed} cells ({removed / initial_count * 100:.1f}%). "
        f"Remaining: {filtered_adata.n_obs} cells"
    )

    return filtered_adata

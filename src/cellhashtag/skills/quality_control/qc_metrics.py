import re
from typing import Any

import numpy as np

try:
    import scanpy as sc
    from anndata import AnnData
except ImportError:
    raise ImportError("scanpy and anndata are required for quality control functions")


def compute_qc_metrics(adata: AnnData) -> dict[str, Any]:
    mito_genes = [
        gene for gene in adata.var_names
        if re.match(r"^MT-|^mt-", gene)
    ]

    ribo_genes = [
        gene for gene in adata.var_names
        if re.match(r"^RPS|^RPL|^RPSA", gene)
    ]

    mito_mask = adata.var_names.isin(mito_genes)
    ribo_mask = adata.var_names.isin(ribo_genes)

    adata.var["mito"] = mito_mask
    adata.var["ribo"] = ribo_mask

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mito", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    adata.obs.rename(columns={
        "n_genes_by_counts": "n_genes",
        "total_counts": "total_counts",
        "pct_counts_mito": "pct_mito",
        "pct_counts_ribo": "pct_ribo"
    }, inplace=True)

    return {
        "n_genes": adata.obs["n_genes"].to_dict(),
        "total_counts": adata.obs["total_counts"].to_dict(),
        "pct_mito": adata.obs["pct_mito"].to_dict(),
        "pct_ribo": adata.obs["pct_ribo"].to_dict()
    }

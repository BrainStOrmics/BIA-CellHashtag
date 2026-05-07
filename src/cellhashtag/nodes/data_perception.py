"""
Data Perception 节点。

感知输入数据的基本信息：物种、平台、细胞数、基因数等。
为后续 Metadata 提取和组学类型判断提供基础。
"""

import pickle
from typing import Any


def node_data_perception(state: dict) -> dict:
    """
    感知输入数据，提取基本信息。

    Args:
        state: CellHashtagState。

    Returns:
        更新后的 state 片段。
    """
    import scanpy as sc

    input_path = state["input_path"]
    input_format = state.get("input_format", "h5ad")

    # 加载数据
    adata = _load_data(input_path, input_format)

    # 感知基本信息
    perceived = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "gene_names": list(adata.var_names[:10]),  # 前 10 个基因名
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "has_raw": adata.raw is not None,
        "is_sparse": hasattr(adata.X, "toarray"),
    }

    # 推断组学类型
    omics_type = _infer_omics_type(adata, state.get("omics_type"))
    perceived["inferred_omics_type"] = omics_type

    # 保存 AnnData 到临时文件
    import tempfile
    adata_path = state.get("adata_path", "adata_temp.pkl")
    with open(adata_path, "wb") as f:
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)

    return {
        "perceived_info": perceived,
        "adata_path": adata_path,
        "omics_type": omics_type,
        "status": "perception_done",
    }


def _load_data(path: str, fmt: str) -> Any:
    """根据格式加载数据。"""
    from ..skills.format_conversion import load_h5ad, load_rds, load_mtx

    if fmt == "h5ad":
        return load_h5ad(path)
    elif fmt == "rds":
        return load_rds(path)
    elif fmt == "mtx":
        return load_mtx(path)
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def _infer_omics_type(adata: Any, provided: str = None) -> str:
    """
    根据 AnnData 内容推断组学类型。

    如果用户已提供，则使用用户提供的值。
    """
    if provided and provided != "auto":
        return provided

    obs_cols = set(adata.obs.columns)
    obsm_keys = set(adata.obsm.keys())

    # 简单的启发式推断
    if "spatial" in obsm_keys or "spatial_coords" in obs_cols:
        return "spRNA"
    if "protein" in obsm_keys or any("adt" in c.lower() for c in obs_cols):
        return "CITE"
    if any("peak" in c.lower() for c in obs_cols) or any("atac" in c.lower() for c in obs_cols):
        return "ATAC"

    return "scRNA"

"""
RDS (Seurat) 格式加载器。

通过 anndata2ri 或 zellkonverter 将 Seurat 对象转换为 AnnData。
需要 R 环境和 rpy2。
"""

import os
from typing import Optional


def load_rds(
    path: str,
    method: str = "anndata2ri",
    seurat_object_name: Optional[str] = None,
) -> any:
    """
    加载 .rds Seurat 对象并转换为 AnnData。

    Args:
        path: .rds 文件路径。
        method: 转换方法，"anndata2ri" 或 "zellkonverter"。
        seurat_object_name: Seurat 对象名称（可选）。

    Returns:
        AnnData 对象。

    Raises:
        ImportError: 如果依赖包未安装。
        FileNotFoundError: 如果文件不存在。
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"RDS file not found: {path}")

    if method == "zellkonverter":
        return _load_via_zellkonverter(path)
    else:
        return _load_via_anndata2ri(path, seurat_object_name)


def _load_via_anndata2ri(path: str, obj_name: Optional[str] = None) -> any:
    """通过 anndata2ri 加载。"""
    try:
        import rpy2.robjects as ro
        import anndata2ri
        anndata2ri.activate()
    except ImportError:
        raise ImportError(
            "rpy2 and anndata2ri are required. "
            "Install with: pip install rpy2 anndata2ri"
        )

    obj_var = obj_name or "seurat_obj"
    ro.r(f'{obj_var} <- readRDS("{path}")')
    ro.r("library(SeuratDisk)")
    ro.r(f"{obj_var} <- as.SingleCellExperiment({obj_var})")

    adata = ro.conversion.rpy2py(ro.globalenv[obj_var])
    return adata


def _load_via_zellkonverter(path: str) -> any:
    """通过 zellkonverter 加载。"""
    try:
        import rpy2.robjects as ro
    except ImportError:
        raise ImportError(
            "rpy2 is required. Install with: pip install rpy2"
        )

    ro.r("library(zellkonverter)")
    ro.r(f'adata <- readH5Seurat("{path}", mode="r")')
    adata = ro.conversion.rpy2py(ro.globalenv["adata"])
    return adata

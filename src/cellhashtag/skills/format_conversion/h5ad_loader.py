"""
H5AD 格式加载器。

AnnData 原生格式，直接加载。
"""

import os
from typing import Any


def load_h5ad(path: str) -> Any:
    """
    加载 .h5ad 文件。

    Args:
        path: .h5ad 文件路径。

    Returns:
        AnnData 对象。

    Raises:
        FileNotFoundError: 如果文件不存在。
        ValueError: 如果文件格式无效。
    """
    import scanpy as sc

    if not os.path.exists(path):
        raise FileNotFoundError(f"H5AD file not found: {path}")

    adata = sc.read_h5ad(path)
    return adata

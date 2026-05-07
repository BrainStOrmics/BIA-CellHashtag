"""
10x MTX 格式加载器。

支持 10x Genomics 的 matrix.mtx + barcodes.tsv + features.tsv 格式。
"""

import os
from typing import Optional


def load_mtx(
    input_dir: str,
    genome: Optional[str] = None,
    gex_only: bool = True,
) -> any:
    """
    加载 10x MTX 格式数据。

    Args:
        input_dir: 包含 matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz 的目录。
        genome: 基因组名称（如 "GRCh38"）。
        gex_only: 仅加载基因表达数据。

    Returns:
        AnnData 对象。

    Raises:
        FileNotFoundError: 如果目录不存在。
    """
    import scanpy as sc

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"MTX directory not found: {input_dir}")

    if genome:
        adata = sc.read_10x_mtx(input_dir, genome=genome, gex_only=gex_only)
    else:
        adata = sc.read_10x_mtx(input_dir, gex_only=gex_only)

    return adata

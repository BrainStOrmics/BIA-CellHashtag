"""SingleCellExperiment .rds loader."""
from typing import Any


def load_sce(path: str, **kwargs: Any) -> Any:
    """Load Seurat/SCE .rds into AnnData.

    Requires rpy2 and zellkonverter. Falls back gracefully if unavailable.

    Args:
        path: Path to .rds file
        **kwargs: Passed to zellkonverter.readRDS

    Returns:
        AnnData object

    Raises:
        ImportError: If rpy2 or zellkonverter not installed
        FileNotFoundError: If path does not exist
    """
    try:
        import zellkonverter
    except ImportError as exc:
        raise ImportError(
            "zellkonverter required for .rds loading. Install: pip install zellkonverter"
        ) from exc

    import os
    if not os.path.exists(path):
        raise FileNotFoundError(f"SCE file not found: {path}")

    return zellkonverter.readRDS(path, **kwargs)

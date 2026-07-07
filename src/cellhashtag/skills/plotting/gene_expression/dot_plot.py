"""Dot plots for gene expression."""
from pathlib import Path
from typing import Any, Optional


def dot_plot(
    adata: Any,
    genes: list[str],
    groupby: str,
    output_dir: Optional[Path] = None,
) -> Optional[str]:
    import scanpy as sc
    import matplotlib.pyplot as plt

    existing = [g for g in genes if g in adata.var_names]
    if not existing:
        raise ValueError(f"No genes found in adata.var_names: {genes}")

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not in adata.obs")

    fig = sc.pl.dotplot(adata, existing, groupby=groupby, swap_axes=True, show=False)

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / f"dotplot_{groupby}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return str(out)

    plt.show()
    return None

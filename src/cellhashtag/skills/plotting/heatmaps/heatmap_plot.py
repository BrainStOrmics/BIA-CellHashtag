"""Heatmap plots for marker genes."""
from pathlib import Path
from typing import Any, Optional


def heatmap(
    adata: Any,
    genes: list[str],
    groupby: str,
    output_dir: Optional[Path] = None,
    dendrogram: bool = True,
    show_values: bool = False,
) -> Optional[str]:
    import scanpy as sc
    import matplotlib.pyplot as plt

    existing = [g for g in genes if g in adata.var_names]
    if not existing:
        raise ValueError(f"No genes found in adata.var_names: {genes}")

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not in adata.obs")

    sc.tl.dendrogram(adata, groupby=groupby, use_rep="X_pca") if dendrogram else None

    fig = sc.pl.heatmap(
        adata, existing, groupby=groupby,
        show_gene_labels=True, dendrogram=dendrogram,
        show=False,
    )

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / f"heatmap_{groupby}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return str(out)

    plt.show()
    return None

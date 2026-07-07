"""PAGA trajectory inference."""
from pathlib import Path
from typing import Any, Optional


def paga_trajectory(
    adata: Any,
    cluster_key: str = "leiden",
    plot: bool = True,
    output_dir: Optional[Path] = None,
) -> dict:
    import scanpy as sc
    import matplotlib.pyplot as plt

    if cluster_key not in adata.obs.columns:
        raise ValueError(f"cluster_key '{cluster_key}' not in adata.obs")

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)

    sc.tl.paga(adata, groups=cluster_key)

    result: dict = {"connectivities_tree": adata.uns["paga"]["connectivities_tree"]}

    if plot:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.paga(adata, ax=ax, show=False)

        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            out = output_dir / "paga_trajectory.png"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            result["plot_path"] = str(out)
        else:
            plt.show()

    return result

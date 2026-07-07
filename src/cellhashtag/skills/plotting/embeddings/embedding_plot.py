"""Publication-quality embedding plots."""
from pathlib import Path
from typing import Any, Optional

OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]


def plot_embedding(
    adata: Any,
    embedding_key: str,
    color_by: Optional[str] = None,
    output_dir: Optional[Path] = None,
    title: Optional[str] = None,
    colors: Optional[list[str]] = None,
) -> Optional[str]:
    import matplotlib.pyplot as plt

    if embedding_key not in adata.obsm:
        raise ValueError(f"Embedding {embedding_key} not found in adata.obsm")

    coords = adata.obsm[embedding_key]
    palette = colors or OKABE_ITO

    fig, ax = plt.subplots(figsize=(8, 6))

    if color_by and color_by in adata.obs.columns:
        values = adata.obs[color_by]
        if values.dtype == object or hasattr(values, "cat"):
            cats = values.astype("category").cat.categories
            cmap = plt.cm.colors.ListedColormap(palette[: len(cats)])
            codes = values.astype("category").cat.codes
            ax.scatter(coords[:, 0], coords[:, 1], c=codes, cmap=cmap, s=5, alpha=0.8)
            for i, cat in enumerate(cats):
                ax.scatter([], [], c=palette[i % len(palette)], label=str(cat), s=30)
            ax.legend(loc="best", markerscale=2, fontsize=6)
        else:
            sc = ax.scatter(coords[:, 0], coords[:, 1], c=values.astype(float),
                            cmap="viridis", s=5, alpha=0.8)
            plt.colorbar(sc, ax=ax, label=color_by)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], c="#999999", s=5, alpha=0.8)

    ax.set_title(title or f"Embedding: {embedding_key}")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    plt.tight_layout()

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        base = f"embedding_{embedding_key}"
        fig.savefig(output_dir / f"{base}.png", dpi=300)
        fig.savefig(output_dir / f"{base}.pdf")
        plt.close(fig)
        return str(output_dir / f"{base}.png")

    plt.show()
    return None

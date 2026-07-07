"""CellHashtag skills library."""

from .clustering_quality import silhouette_check, modularity_check
from .format_conversion import load_h5ad, load_rds, load_mtx, load_sce
from .quality_control import compute_qc_metrics, filter_cells, detect_doublets

__all__ = [
    "silhouette_check",
    "modularity_check",
    "load_h5ad",
    "load_rds",
    "load_mtx",
    "load_sce",
    "compute_qc_metrics",
    "filter_cells",
    "detect_doublets",
]

try:
    from .clustering_quality import phiclust_check
    __all__.append("phiclust_check")
except ImportError:
    phiclust_check = None

try:
    from .plotting.embeddings import plot_embedding
    from .plotting.gene_expression import violin_plot, dot_plot
    from .plotting.heatmaps import heatmap
    __all__.extend(["plot_embedding", "violin_plot", "dot_plot", "heatmap"])
except ImportError:
    pass

try:
    from .clustering_degs import find_markers, consensus_markers, resolution_scan
    __all__.extend(["find_markers", "consensus_markers", "resolution_scan"])
except ImportError:
    pass

try:
    from .batch_correction import correct_harmony, correct_bbknn, correct_combat
    __all__.extend(["correct_harmony", "correct_bbknn", "correct_combat"])
except ImportError:
    pass

try:
    from .spatial import load_visium, load_merfish, load_stereoseq, load_spatial
    __all__.extend(["load_visium", "load_merfish", "load_stereoseq", "load_spatial"])
except ImportError:
    pass

try:
    from .trajectory import paga_trajectory
    __all__.append("paga_trajectory")
except ImportError:
    pass

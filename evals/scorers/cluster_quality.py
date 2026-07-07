"""Cluster quality scorers — structural metrics on an already-clustered AnnData."""
from __future__ import annotations


def _bucket(value: float, thresholds: list[tuple[float, float]]) -> float:
    for lo, score in thresholds:
        if value >= lo:
            return score
    return thresholds[-1][1]


def silhouette_score(fixture, actual: dict) -> float:
    asw = float(actual.get("asw", 0.0))
    n_clusters = int(actual.get("n_clusters", 0))
    if asw > 0.70 and n_clusters > 20:
        return 0.50
    if asw >= 0.50:
        return 1.0
    if asw >= 0.35:
        return 0.85
    if asw >= 0.25:
        return 0.60
    return 0.20


def modularity_score(fixture, actual: dict) -> float:
    mod = float(actual.get("modularity", 0.0))
    if mod >= 0.60:
        return 1.0
    if mod >= 0.40:
        return 0.85
    if mod >= 0.20:
        return 0.55
    return 0.15


def size_distribution_score(fixture, actual: dict) -> float:
    sizes = actual.get("cluster_sizes", [])
    if not sizes:
        return 0.0
    ratio = max(sizes) / max(min(sizes), 1)
    min_cells = fixture.expected.get("min_cells_per_cluster", 20)
    has_tiny = any(s < min_cells for s in sizes)
    if ratio <= 10:
        base = 1.0
    elif ratio <= 30:
        base = 0.75
    elif ratio <= 50:
        base = 0.45
    else:
        base = 0.15
    return base * (0.5 if has_tiny else 1.0)


def marker_coherence_score(fixture, actual: dict) -> float:
    coherence = float(actual.get("lineage_coherence", 0.0))
    if coherence >= 0.80:
        return 1.0
    if coherence >= 0.60:
        return 0.70
    return 0.30

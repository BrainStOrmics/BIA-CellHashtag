"""Annotation scorers — deterministic comparisons against gold cell types."""
from __future__ import annotations


def _normalize(name: str) -> str:
    return " ".join(name.lower().replace("+", " ").replace("-", " ").split())


def exact_match(fixture, actual: dict) -> float:
    predicted = _normalize(actual.get("cell_type", ""))
    gold = _normalize(fixture.expected["cell_type"])
    alts = [_normalize(a) for a in fixture.expected.get("acceptable_alternates", [])]
    excluded = [_normalize(e) for e in fixture.expected.get("excluded_types", [])]

    if predicted in excluded:
        return 0.0
    if predicted == gold or predicted in alts:
        return 1.0
    if actual.get("lineage", "").lower() == fixture.expected.get("lineage", "").lower():
        return 0.50
    if predicted == "unknown":
        return 0.10
    return 0.0


def tier1_jaccard(fixture, actual: dict) -> float:
    gold = set(fixture.expected.get("tier1_markers", []))
    predicted = set(actual.get("evidence_markers", []))
    if not gold:
        return 1.0 if not predicted else 0.0
    inter = len(gold & predicted)
    union = len(gold | predicted)
    wj = inter / union if union else 0.0
    if wj >= 0.60:
        return 1.0
    if wj >= 0.40:
        return 0.75
    if wj >= 0.20:
        return 0.45
    return 0.10


def expected_calibration_error(fixture, actual: dict, n_bins: int = 10) -> float:
    """Returns 1 - ECE so higher is better (ECE=0 → score=1.0)."""
    conf = float(actual.get("confidence", 0.5))
    correct = 1.0 if exact_match(fixture, actual) >= 1.0 else 0.0
    ece = abs(conf - correct)
    if ece <= 0.05:
        return 1.0
    if ece <= 0.10:
        return 0.80
    if ece <= 0.15:
        return 0.55
    return 0.15


def excluded_avoidance(fixture, actual: dict) -> float:
    excluded = {_normalize(e) for e in fixture.expected.get("excluded_types", [])}
    top3 = [_normalize(t) for t in actual.get("top3_predictions", [])]
    hits = sum(1 for t in top3 if t in excluded)
    if hits == 0:
        return 1.0
    if hits == 1:
        return 0.30
    return 0.0

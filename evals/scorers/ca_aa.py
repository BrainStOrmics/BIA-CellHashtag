"""CA-AA alternation decision scorers."""
from __future__ import annotations


def action_match(fixture, actual: dict) -> float:
    exp = fixture.expected["action"]
    got = actual.get("action")
    if got == exp:
        return 1.0
    if {got, exp} == {"merge", "subcluster"}:
        return 0.20
    if got == "converge" and exp != "converge":
        return 0.10
    return 0.0


def _pairs_equal_unordered(a: list, b: list) -> bool:
    norm = lambda xs: sorted(tuple(sorted(p)) for p in xs)
    return norm(a or []) == norm(b or [])


def merge_correctness(fixture, actual: dict) -> float:
    exp = fixture.expected.get("merge_candidates", [])
    got = actual.get("merge_candidates", [])
    if _pairs_equal_unordered(exp, got):
        return 1.0
    if not exp and not got:
        return 1.0
    exp_flat = {tuple(sorted(p)) for p in exp}
    got_flat = {tuple(sorted(p)) for p in got}
    if exp_flat and got_flat:
        overlap = len(exp_flat & got_flat)
        if overlap:
            return 0.50
    return 0.0


def subcluster_correctness(fixture, actual: dict) -> float:
    exp = set(fixture.expected.get("subcluster_targets", []))
    got = set(actual.get("subcluster_targets", []))
    if exp == got:
        return 1.0
    if exp & got:
        return 0.50
    if exp and not got:
        return 0.20
    return 0.0


def convergence_signal(fixture, actual: dict) -> float:
    exp = fixture.expected.get("converged", False)
    got = actual.get("converged", False)
    if exp == got:
        return 1.0
    if got and not exp:
        return 0.0
    return 0.30

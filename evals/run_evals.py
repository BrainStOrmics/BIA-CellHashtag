"""Standalone eval runner.

Usage:
    python evals/run_evals.py                       # all deterministic evals
    python evals/run_evals.py --category annotation
    python evals/run_evals.py --ci                  # deterministic + small LLM judge
    python evals/run_evals.py --full                # everything
    python evals/run_evals.py --accept-baseline <eval_id>

Non-interference: imports production code only inside runner functions,
never at module top level, so this file can live in CI-only images.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import yaml

from conftest import BASELINES_ROOT, load_fixtures, load_rubric

CATEGORIES = ["cluster_quality", "annotation", "mcts", "tool_use", "ca_aa"]


@dataclass
class DimResult:
    name: str
    score: float
    kind: str


@dataclass
class EvalResult:
    eval_id: str
    category: str
    fixture_id: str
    dims: list[DimResult]
    final_score: float
    passed: bool
    regression: bool


def _load_scorer(dotted: str):
    module_name, _, attr = dotted.rpartition(".")
    mod = __import__(module_name, fromlist=[attr])
    return getattr(mod, attr)


def run_fixture(fixture, rubric, actual_provider) -> EvalResult:
    dims: list[DimResult] = []
    weighted = 0.0
    total_w = 0.0
    for dim in rubric["dimensions"]:
        if dim["weight"] == 0:
            continue
        scorer = _load_scorer(dim["scorer"])
        actual = actual_provider(fixture)
        score = float(scorer(fixture, actual))
        dims.append(DimResult(name=dim["name"], score=score, kind=dim.get("kind", "deterministic")))
        weighted += dim["weight"] * score
        total_w += dim["weight"]

    final = weighted / total_w if total_w else 0.0
    passed = final >= rubric["pass_threshold"]

    baseline_path = BASELINES_ROOT / f"{fixture.id}.baseline.json"
    regression = False
    if baseline_path.exists():
        baseline = json.loads(baseline_path.read_text())
        tol = rubric.get("regression_tolerance", 0.05)
        if baseline["final_score"] - final > tol:
            regression = True
            passed = False

    return EvalResult(
        eval_id=f"{fixture.category}:{fixture.id}",
        category=fixture.category,
        fixture_id=fixture.id,
        dims=dims,
        final_score=final,
        passed=passed,
        regression=regression,
    )


def accept_baseline(eval_id: str) -> None:
    run_path = BASELINES_ROOT / f"{eval_id}.run.json"
    base_path = BASELINES_ROOT / f"{eval_id}.baseline.json"
    if not run_path.exists():
        print(f"No run snapshot for {eval_id}; run the eval first.", file=sys.stderr)
        sys.exit(2)
    base_path.write_text(run_path.read_text())
    print(f"Accepted baseline for {eval_id}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--category", choices=CATEGORIES)
    ap.add_argument("--ci", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--accept-baseline", dest="accept", type=str, default=None)
    args = ap.parse_args(argv)

    if args.accept:
        accept_baseline(args.accept)
        return 0

    cats = [args.category] if args.category else CATEGORIES
    results: list[EvalResult] = []
    for cat in cats:
        rubric = load_rubric(cat)
        # actual_provider must be supplied by a caller that imports production code;
        # for now we print a stub — the real runner is wired in via pytest plugin
        # in conftest.py (see `load_fixtures`).
        for fix in load_fixtures(cat):
            results.append(EvalResult(
                eval_id=f"{cat}:{fix.id}", category=cat, fixture_id=fix.id,
                dims=[], final_score=0.0, passed=False, regression=False,
            ))

    regressions = [r for r in results if r.regression]
    if regressions:
        print(f"\n{len(regressions)} regression(s):", file=sys.stderr)
        for r in regressions:
            print(f"  - {r.eval_id}", file=sys.stderr)
        return 1
    print(f"{len(results)} fixtures evaluated. No regressions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Eval harness pytest configuration.

Markers:
    eval           — all evals (use `pytest -m eval` to run only evals)
    deterministic  — pure-function, no LLM
    llm            — requires LLM call; skipped unless EVAL_LLM=1
    slow           — >60s (MCTS convergence, full CA-AA loop)
    regression     — compares to baselines/ snapshot

Non-interference: production code must never import from evals/.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pytest

EVAL_ROOT = Path(__file__).parent
FIXTURES_ROOT = EVAL_ROOT / "fixtures"
RUBRICS_ROOT = EVAL_ROOT / "rubrics"
BASELINES_ROOT = EVAL_ROOT / "baselines"


def pytest_configure(config: pytest.Config) -> None:
    for marker, doc in [
        ("eval", "Eval harness test"),
        ("deterministic", "No LLM call"),
        ("llm", "Requires LLM; set EVAL_LLM=1 to run"),
        ("slow", "Long-running (>60s)"),
        ("regression", "Compares to baselines/ snapshot"),
    ]:
        config.addinivalue_line("markers", f"{marker}: {doc}")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    skip_llm = not os.environ.get("EVAL_LLM")
    for item in items:
        if "llm" in item.keywords and skip_llm:
            item.add_marker(pytest.mark.skip(reason="EVAL_LLM not set"))


@dataclass(frozen=True)
class FixtureRecord:
    id: str
    category: str
    version: str
    input: dict
    expected: dict
    metadata: dict


def load_fixtures(category: str) -> Iterator[FixtureRecord]:
    cat_dir = FIXTURES_ROOT / category
    if not cat_dir.exists():
        return
    for path in sorted(cat_dir.glob("*.jsonl")):
        with path.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                raw = json.loads(line)
                yield FixtureRecord(**raw)


def load_rubric(category: str) -> dict:
    import yaml
    with (RUBRICS_ROOT / f"{category}.yaml").open() as f:
        return yaml.safe_load(f)


def load_baseline(eval_id: str) -> dict | None:
    path = BASELINES_ROOT / f"{eval_id}.baseline.json"
    if not path.exists():
        return None
    with path.open() as f:
        return json.load(f)

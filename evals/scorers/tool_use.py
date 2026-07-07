"""Tool-use scorers — trace-based, deterministic."""
from __future__ import annotations


def correct_tool(fixture, actual: dict) -> float:
    return 1.0 if actual.get("tool_called") == fixture.expected["tool"] else 0.0


def required_args(fixture, actual: dict) -> float:
    expected = fixture.expected.get("required_args", {})
    got = actual.get("args", {})
    for k, v in expected.items():
        if got.get(k) != v:
            return 0.0
    return 1.0


def forbidden_calls(fixture, actual: dict) -> float:
    must_not = set(fixture.expected.get("must_not_call", []))
    called = set(actual.get("trace", []))
    return 0.0 if must_not & called else 1.0


def coverage(fixture, actual: dict) -> float:
    expected_set = set(fixture.expected.get("tools_called_set", []))
    called = set(actual.get("trace", []))
    min_tools = fixture.expected.get("min_tools_called", len(expected_set))
    overlap = expected_set & called
    if len(overlap) >= min_tools:
        return 1.0
    return len(overlap) / max(min_tools, 1)

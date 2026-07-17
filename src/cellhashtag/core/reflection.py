"""LATS structured reflection schema.

LLM outputs structured diagnosis text at MCTS evaluation phase. Fields feed
back into reward computation, prior biasing, and confidence discounting.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from typing import Any, Optional


@dataclass
class ReflectionSchema:
    diagnosis: str = ""
    assessment: float = 0.0
    suggested_next: str = ""
    hallucination_risk: str = "low"
    should_continue: bool = True
    raw: dict = field(default_factory=dict)

    def risk_weight(self) -> float:
        return {"low": 1.0, "medium": 0.6, "high": 0.2}.get(self.hallucination_risk, 0.5)


def parse_reflection(llm_output: Any) -> Optional[ReflectionSchema]:
    content = _extract_text(llm_output)
    if not content:
        return None
    m = re.search(r"\{.*\}", content, re.DOTALL)
    if not m:
        return None
    try:
        data = json.loads(m.group())
    except json.JSONDecodeError:
        return None

    risk = str(data.get("hallucination_risk", "low")).lower()
    if risk not in {"low", "medium", "high"}:
        risk = "medium"

    try:
        assessment = float(data.get("assessment", data.get("action_quality", 0.0)))
    except (TypeError, ValueError):
        assessment = 0.0

    return ReflectionSchema(
        diagnosis=str(data.get("diagnosis", "")),
        assessment=assessment,
        suggested_next=str(data.get("suggested_next", "")),
        hallucination_risk=risk,
        should_continue=bool(data.get("should_continue", True)),
        raw=data,
    )


def reflection_to_dict(r: ReflectionSchema) -> dict:
    return asdict(r)


def _extract_text(result: Any) -> str:
    if isinstance(result, str):
        return result
    messages = result.get("messages", []) if isinstance(result, dict) else []
    for msg in reversed(messages):
        content = getattr(msg, "content", "") if hasattr(msg, "content") else str(msg)
        if content:
            return content
    return ""

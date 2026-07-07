"""Adversarial LLM judge for annotation reasoning.

Used only by the `adversarial_reasoning_quality` dimension in rubrics/annotation.yaml.
MUST:
  - use a different model or temperature than the system under test
  - prompt be adversarial ("find why this is wrong")
  - be calibrated against fixtures/annotation/calibration.jsonl
  - be gated on deterministic dims (see rubric `gated_by`)

Non-interference: imports only stdlib + the project LLM client.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

ADVERSARIAL_PROMPT = """You are auditing a single-cell annotation decision.
Find reasons the annotation could be WRONG.

Cluster markers (top-15): {markers}
Tissue context: {tissue}
Predicted cell type: {predicted}
Agent's reasoning: {reasoning}

Critique in three parts:
1. Markers that contradict the prediction (specific gene names)
2. Alternative cell types that also fit these markers
3. Tissue-context implausibilities

Final score 0.0-1.0: 1.0 = reasoning holds up; 0.0 = reasoning is fatally flawed.
Return ONLY a JSON object: {{"score": <float>, "critique": "<one line>"}}"""


@dataclass(frozen=True)
class JudgeResult:
    score: float
    critique: str
    calibrated: bool


def reasoning_score(fixture, actual: dict, *, llm_client) -> JudgeResult:
    if not os.environ.get("EVAL_LLM"):
        return JudgeResult(score=0.0, critique="LLM disabled", calibrated=False)

    prompt = ADVERSARIAL_PROMPT.format(
        markers=", ".join(fixture.input.get("top_markers", [])),
        tissue=fixture.input.get("tissue", "unknown"),
        predicted=actual.get("cell_type", "Unknown"),
        reasoning=actual.get("reasoning", ""),
    )
    response = llm_client.complete(prompt, temperature=0.0)
    score = float(response["score"])
    return JudgeResult(score=score, critique=response["critique"], calibrated=True)

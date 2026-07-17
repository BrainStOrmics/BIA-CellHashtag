"""Five-component multi-objective reward + sliding-window normalization + potential shaping.

R_raw = w1*R_purity + w2*R_specificity + w3*R_context + w4*R_lats + w5*R_known_marker - Omega
R_tilde = (R_raw - mu_W) / (sigma_W + eps)
Phi(s) = sum_l [w1*Sil(l) + w2*MarkerCoherence(l)] - lambda_frag*|L|
R_final = R_tilde + gamma*Phi(s') - Phi(s)
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass
from typing import Optional


DEFAULT_REWARD_WEIGHTS = {
    "purity": 0.25,
    "specificity": 0.20,
    "context": 0.20,
    "lats": 0.20,
    "known_marker": 0.15,
}
FRAG_PENALTY = 0.02
WINDOW_SIZE = 50
EPS = 1e-8


@dataclass
class RewardBreakdown:
    purity: float = 0.0
    specificity: float = 0.0
    context: float = 0.0
    lats: float = 0.0
    known_marker: float = 0.0
    omega: float = 0.0
    raw: float = 0.0
    normalized: float = 0.0
    potential_delta: float = 0.0
    final: float = 0.0


def compute_reward(
    *,
    purity: float,
    specificity: float,
    context_score: float,
    lats_score: float,
    known_marker_score: float,
    overfragmentation_penalty: float = 0.0,
    weights: Optional[dict[str, float]] = None,
) -> RewardBreakdown:
    w = weights or DEFAULT_REWARD_WEIGHTS
    raw = (
        w.get("purity", 0.0) * purity
        + w.get("specificity", 0.0) * specificity
        + w.get("context", 0.0) * context_score
        + w.get("lats", 0.0) * lats_score
        + w.get("known_marker", 0.0) * known_marker_score
        - overfragmentation_penalty
    )
    return RewardBreakdown(
        purity=purity,
        specificity=specificity,
        context=context_score,
        lats=lats_score,
        known_marker=known_marker_score,
        omega=overfragmentation_penalty,
        raw=raw,
    )


class SlidingWindowNormalizer:
    def __init__(self, window_size: int = WINDOW_SIZE):
        self._buf: deque[float] = deque(maxlen=window_size)

    def push(self, raw_reward: float) -> None:
        self._buf.append(raw_reward)

    def normalize(self, raw_reward: float) -> float:
        if len(self._buf) < 2:
            return raw_reward
        mean = sum(self._buf) / len(self._buf)
        var = sum((x - mean) ** 2 for x in self._buf) / len(self._buf)
        std = math.sqrt(var) if var > 0 else EPS
        return (raw_reward - mean) / (std + EPS)

    def stats(self) -> dict:
        if not self._buf:
            return {"mu": 0.0, "sigma": 0.0, "n": 0}
        mean = sum(self._buf) / len(self._buf)
        var = sum((x - mean) ** 2 for x in self._buf) / len(self._buf)
        return {"mu": mean, "sigma": math.sqrt(var), "n": len(self._buf)}


@dataclass
class ClusterStats:
    label: str
    silhouette: float = 0.0
    marker_coherence: float = 0.0


def potential(state_clusters: list[ClusterStats], lambda_frag: float = FRAG_PENALTY) -> float:
    phi = 0.0
    for c in state_clusters:
        phi += 0.5 * c.silhouette + 0.5 * c.marker_coherence
    phi -= lambda_frag * len(state_clusters)
    return phi


def final_reward(
    breakdown: RewardBreakdown,
    normalized_raw: float,
    phi_next: float,
    phi_current: float,
    gamma: float = 0.95,
) -> RewardBreakdown:
    breakdown.normalized = normalized_raw
    breakdown.potential_delta = gamma * phi_next - phi_current
    breakdown.final = normalized_raw + breakdown.potential_delta
    return breakdown


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa or not sb:
        return 0.0
    inter = len(sa & sb)
    union = len(sa | sb)
    return inter / union if union else 0.0

"""Three-stage cold start: metadata -> profiling -> LLM prior injection.

Builds the root DAGNode before MCTS blind search begins. Shrinks initial
search space by seeding from tissue templates + quick scanpy cluster + LLM priors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional


TISSUE_TEMPLATES: dict[str, dict[str, Any]] = {
    "pbmc": {"expected_clusters": (8, 12), "resolution": 0.5, "lineages": ["T", "B", "NK", "Mono", "DC", "Platelet"]},
    "tumor": {"expected_clusters": (15, 25), "resolution": 0.8, "lineages": ["Tumor", "T", "B", "Macro", "Fibro", "Endo"]},
    "brain": {"expected_clusters": (10, 20), "resolution": 0.6, "lineages": ["Neuron", "Astro", "Micro", "Oligo", "Endo"]},
    "lung": {"expected_clusters": (10, 18), "resolution": 0.6, "lineages": ["Epithelial", "T", "Macro", "Fibro", "Endo"]},
    "liver": {"expected_clusters": (8, 14), "resolution": 0.5, "lineages": ["Hepato", "Kupffer", "Endo", "T", "NK"]},
    "default": {"expected_clusters": (8, 20), "resolution": 0.5, "lineages": []},
}


@dataclass
class ColdStartResult:
    tissue: str
    resolution_hint: float
    expected_lineages: list[str]
    quick_clusters: list[str]
    llm_priors: list[dict] = field(default_factory=list)
    action_shortlist: list[dict] = field(default_factory=list)


def match_tissue(metadata: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    tissue = str(metadata.get("tissue", "")).lower()
    for key in TISSUE_TEMPLATES:
        if key != "default" and key in tissue:
            return key, TISSUE_TEMPLATES[key]
    return "default", TISSUE_TEMPLATES["default"]


def quick_profile(adata: Any, resolution: float, cluster_key: str = "leiden_quick") -> list[str]:
    try:
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(1500, adata.n_vars))
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)
        return sorted(set(adata.obs[cluster_key].astype(str).tolist()))
    except Exception:
        return []


def llm_prior_prompt(tissue: str, lineages: list[str], quick_clusters: list[str], markers_per_cluster: dict) -> str:
    return (
        f"Tissue: {tissue}\n"
        f"Expected lineages: {', '.join(lineages) or 'unknown'}\n"
        f"Quick clusters: {', '.join(quick_clusters)}\n"
        f"Per-cluster top markers:\n"
        + "\n".join(f"  {k}: {', '.join(v[:10])}" for k, v in markers_per_cluster.items())
        + "\nFor each cluster, propose initial cell type label + confidence + action suggestion "
          "(merge/split/freeze). Output JSON: "
          '{"priors": [{"cluster": "0", "label": "...", "confidence": 0.0-1.0, "action": "merge|split|freeze|none", "target": "..."}]}'
    )


def cold_start(
    metadata: dict[str, Any],
    adata: Optional[Any] = None,
    llm_fn: Optional[Any] = None,
    markers_per_cluster: Optional[dict] = None,
) -> ColdStartResult:
    tissue_name, template = match_tissue(metadata)
    resolution = template["resolution"]
    lineages = list(template["lineages"])

    quick_clusters: list[str] = []
    if adata is not None:
        quick_clusters = quick_profile(adata, resolution)

    priors: list[dict] = []
    shortlist: list[dict] = []
    if llm_fn is not None and quick_clusters and markers_per_cluster:
        try:
            prompt = llm_prior_prompt(tissue_name, lineages, quick_clusters, markers_per_cluster)
            raw = llm_fn(prompt)
            priors, shortlist = _parse_llm_priors(raw)
        except Exception:
            priors, shortlist = [], []

    return ColdStartResult(
        tissue=tissue_name,
        resolution_hint=resolution,
        expected_lineages=lineages,
        quick_clusters=quick_clusters,
        llm_priors=priors,
        action_shortlist=shortlist,
    )


def _parse_llm_priors(raw: Any) -> tuple[list[dict], list[dict]]:
    import json
    import re
    text = raw if isinstance(raw, str) else str(getattr(raw, "content", raw))
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if not m:
        return [], []
    try:
        data = json.loads(m.group())
    except json.JSONDecodeError:
        return [], []
    priors = data.get("priors", []) if isinstance(data, dict) else []
    shortlist = [
        {"action": p.get("action"), "cluster": p.get("cluster"), "target": p.get("target"), "weight": p.get("confidence", 0.5)}
        for p in priors if p.get("action") in {"merge", "split", "freeze"}
    ]
    return priors, shortlist

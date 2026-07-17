"""CellHashtag configuration with Pydantic validation, profile presets, and env override."""

import os
from pathlib import Path
from typing import Any, Optional

try:
    import yaml
except ImportError:
    yaml = None

from pydantic import BaseModel, Field


# ===== Profile presets =====

PROFILES: dict[str, dict[str, Any]] = {
    "fast": {
        "llm": {"temperature": 0.1, "max_tokens": 1024},
        "clustering": {"max_iterations": 1},
        "annotation": {"max_anno_iter": 2},
        "lats": {
            "search_params": {"n_iterations": 3, "max_branches": 1},
        },
    },
    "default": {},
    "deep": {
        "llm": {"temperature": 0.2, "max_tokens": 4096},
        "clustering": {"max_iterations": 5},
        "annotation": {"max_anno_iter": 8},
        "lats": {
            "search_params": {"n_iterations": 20, "max_branches": 5},
        },
    },
}

# ===== Environment variable mappings =====
# Keys are dotted paths into the config dict, values are env var names.
_ENV_OVERRIDES: dict[str, str] = {
    "llm.api_key": "CELLHASHTAG_LLM_API_KEY",
    "llm.provider": "CELLHASHTAG_LLM_PROVIDER",
    "llm.model": "CELLHASHTAG_LLM_MODEL",
    "llm.api_base": "CELLHASHTAG_LLM_API_BASE",
    "llm.temperature": "CELLHASHTAG_LLM_TEMPERATURE",
    "llm.max_tokens": "CELLHASHTAG_LLM_MAX_TOKENS",
    "llm.timeout": "CELLHASHTAG_LLM_TIMEOUT",
    "llm.enable_thinking": "CELLHASHTAG_LLM_ENABLE_THINKING",
    "llm.enable_search": "CELLHASHTAG_LLM_ENABLE_SEARCH",
    "clustering.max_iterations": "CELLHASHTAG_CLUSTERING_MAX_ITERATIONS",
    "clustering.default_resolution": "CELLHASHTAG_CLUSTERING_DEFAULT_RESOLUTION",
    "annotation.max_anno_iter": "CELLHASHTAG_ANNOTATION_MAX_ANNO_ITER",
    "annotation.confidence_threshold": "CELLHASHTAG_ANNOTATION_CONFIDENCE_THRESHOLD",
    "annotation.marker_top_k": "CELLHASHTAG_ANNOTATION_MARKER_TOP_K",
    "lats.search_params.n_iterations": "CELLHASHTAG_LATS_N_ITERATIONS",
    "lats.search_params.max_branches": "CELLHASHTAG_LATS_MAX_BRANCHES",
    "lats.search_params.confidence_threshold": "CELLHASHTAG_LATS_CONFIDENCE_THRESHOLD",
    "cellwiki.wiki_dir": "CELLHASHTAG_CELLWIKI_DIR",
    "output.dir": "CELLHASHTAG_OUTPUT_DIR",
}


def _set_nested(data: dict, dotted_path: str, value: Any):
    keys = dotted_path.split(".")
    current = data
    for k in keys[:-1]:
        if k not in current:
            current[k] = {}
        current = current[k]
    current[keys[-1]] = value


def _deep_merge(base: dict, override: dict) -> dict:
    result = base.copy()
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = val
    return result


def _coerce_env_value(raw: str, target_type: type | None = None) -> Any:
    if target_type is bool:
        return raw.lower() in ("true", "1", "yes", "on")
    if target_type is int:
        return int(raw)
    if target_type is float:
        return float(raw)
    return raw


def _get_type_hints(model_cls, path: str) -> type | None:
    try:
        keys = path.split(".")
        cls = model_cls
        for k in keys:
            ann = cls.__annotations__.get(k)
            if hasattr(ann, "__origin__"):
                cls = ann.__origin__
            elif isinstance(ann, type):
                cls = ann
            else:
                return None
        return cls
    except Exception:
        return None


def _collect_env_overrides() -> dict:
    overrides: dict = {}
    for path, env_var in _ENV_OVERRIDES.items():
        val = os.environ.get(env_var)
        if val is not None:
            hint = _get_type_hints(CellHashtagConfig, path)
            overrides[path] = _coerce_env_value(val, hint)
    return overrides


def _apply_env_overrides(raw: dict) -> dict:
    for path, val in _collect_env_overrides().items():
        _set_nested(raw, path, val)
    return raw


# ===== Pydantic models =====

class LLMConfig(BaseModel):
    provider: str = "dashscope"
    model: str = "qwen3.5-plus"
    api_key: str = ""
    api_base: str = "https://dashscope.aliyuncs.com/compatible-mode/v1"
    temperature: float = 0.3
    max_tokens: int = 2048
    timeout: int = 120
    enable_thinking: bool = True
    enable_search: bool = False

    def get_api_key(self) -> str:
        if self.api_key:
            return self.api_key
        return os.environ.get("CELLHASHTAG_LLM_API_KEY", "")


class ClusteringConfig(BaseModel):
    default_resolution: float = 0.5
    resolution_range: tuple[float, float] = (0.05, 2.0)
    quality_metrics: list[str] = Field(default_factory=lambda: ["silhouette_score", "modularity"])
    min_cluster_size: int = 10
    max_clusters: int = 50
    max_iterations: int = 3


class AnnotationConfig(BaseModel):
    max_anno_iter: int = 5
    confidence_threshold: float = 0.7
    marker_top_k: int = 20


class LATSSearchConfig(BaseModel):
    n_iterations: int = 10
    exploration_weight: float = 1.414
    max_branches: int = 3
    confidence_threshold: float = 0.7
    early_stop_threshold: float = 0.95


class LATSConfig(BaseModel):
    search_params: LATSSearchConfig = Field(default_factory=LATSSearchConfig)
    value_weights: dict[str, float] = Field(default_factory=lambda: {
        "marker_match": 0.4,
        "ontology_consistency": 0.3,
        "evidence_diversity": 0.2,
        "llm_confidence": 0.1,
    })


class DAGMCTSConfig(BaseModel):
    max_iterations: int = 20
    cpuct: float = 1.5
    tau0: float = 1.0
    tau_gamma: float = 0.25
    alpha: float = 0.3
    epsilon: float = 0.05
    early_stop_threshold: float = 0.95
    confidence_threshold: float = 0.7
    fpu_delta: float = 0.5
    theta_high: float = 0.7
    theta_low: float = 0.3
    reward_weights: dict[str, float] = Field(default_factory=lambda: {
        "purity": 0.25,
        "specificity": 0.20,
        "context": 0.20,
        "lats": 0.20,
        "known_marker": 0.15,
    })
    enable_tree_viz: bool = True


class CellWikiConfig(BaseModel):
    wiki_dir: str = "~/CellWiki/wiki"

    def get_wiki_dir(self) -> Path:
        return Path(self.wiki_dir).expanduser()


class OutputConfig(BaseModel):
    dir: str = "output"
    report_filename: str = "annotation_report.md"
    table_filename: str = "annotation_table.csv"
    include_evidence: bool = True


class CellHashtagConfig(BaseModel):
    llm: LLMConfig = Field(default_factory=LLMConfig)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    lats: LATSConfig = Field(default_factory=LATSConfig)
    dag_mcts: DAGMCTSConfig = Field(default_factory=DAGMCTSConfig)
    cellwiki: CellWikiConfig = Field(default_factory=CellWikiConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)


def load_config(
    path: Optional[str] = None,
    *,
    profile: str = "default",
    **overrides,
) -> CellHashtagConfig:
    """Load config with three-layer override.

    Priority (lowest → highest):
      1. Built-in defaults in Pydantic models
      2. YAML file (default: config.yaml next to this module)
      3. Environment variables (CELLHASHTAG_*)
      4. Profile preset (fast / default / deep)
      5. Runtime kwargs

    Args:
        path: YAML file path. None → default path.
        profile: Named preset ("fast", "default", "deep").
        **overrides: Runtime overrides, e.g. llm__model="gpt-4o".

    Returns:
        Validated CellHashtagConfig.

    Raises:
        ValueError: If config validation fails (with clear error message).
    """
    if yaml is None:
        raise ImportError("PyYAML is required: pip install pyyaml")

    if path is None:
        path = str(Path(__file__).parent / "config.yaml")

    raw: dict = {}

    # Layer 1: YAML file
    if os.path.exists(path):
        with open(path) as f:
            raw = yaml.safe_load(f) or {}
    else:
        import warnings
        warnings.warn(
            f"Config file not found: {path}. "
            f"Using defaults. Set CELLHASHTAG_LLM_API_KEY or create config.yaml.",
            UserWarning,
            stacklevel=2,
        )

    # Layer 2: Environment variables
    raw = _apply_env_overrides(raw)

    # Layer 3: Profile preset
    if profile != "default" and profile in PROFILES:
        raw = _deep_merge(raw, PROFILES[profile])

    # Layer 4: Runtime kwargs (use __ to nest, e.g. llm__model="gpt-4o")
    for key, val in overrides.items():
        dotted = key.replace("__", ".")
        _set_nested(raw, dotted, val)

    try:
        return CellHashtagConfig(**raw)
    except Exception as e:
        raise ValueError(f"Invalid config: {e}") from e


def setup_llm(config: LLMConfig):
    """Create a LangChain ChatOpenAI instance from config."""
    try:
        from langchain_openai import ChatOpenAI
    except ImportError:
        print("Warning: langchain-openai not installed.")
        return None

    api_key = config.get_api_key()
    if not api_key:
        raise ValueError(
            "LLM API key not configured. Set CELLHASHTAG_LLM_API_KEY env var "
            "or add 'api_key' to config.yaml."
        )

    return ChatOpenAI(
        api_key=api_key,
        base_url=config.api_base,
        model=config.model,
        temperature=config.temperature,
        max_tokens=config.max_tokens,
        request_timeout=config.timeout,
        extra_body={
            "enable_thinking": config.enable_thinking,
            "enable_search": config.enable_search,
        },
    )

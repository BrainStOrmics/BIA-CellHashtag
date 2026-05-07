"""
CellHashtag 配置加载

保留现有 config.yaml 格式，支持 LLM API 配置和 Agent 参数。
"""

import os
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# ===== LLM 配置 =====
class LLMConfig:
    """LLM API 配置。"""
    CHAT_MODEL_API = {
        "api": "",
        "url": "",
        "model": "",
        "type": "openai",
    }
    MULTIMODAL_CHAT_MODEL_API = {
        "api": "",
        "url": "",
        "model": "",
        "type": "openai",
    }
    ENABLE_THINKING = True
    ENABLE_SEARCH = False


# ===== Agent 配置 =====
class AgentConfig:
    """Agent 运行参数。"""
    DEFAULT_CLUSTER_KEY = "leiden"
    DEFAULT_MAX_ITERATIONS = 3         # 聚类迭代上限
    DEFAULT_MAX_ANNO_ITER = 5          # 注释自批评迭代上限
    DEFAULT_OMICS_TYPE = "scRNA"
    OUTPUT_DIR = "output"
    REPORT_FILENAME = "annotation_report.md"
    TABLE_FILENAME = "annotation_table.csv"


def load_config(yaml_path: str = None) -> dict:
    """
    从 config.yaml 加载配置。

    Args:
        yaml_path: YAML 配置文件路径。如果为 None，使用默认路径。

    Returns:
        配置字典。
    """
    if yaml is None:
        raise ImportError("PyYAML is required. Install with: pip install pyyaml")

    if yaml_path is None:
        yaml_path = str(Path(__file__).parent / "config.yaml")

    if not os.path.exists(yaml_path):
        print(f"Warning: Config file not found at {yaml_path}, using defaults.")
        return {}

    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f) or {}

    # 应用 LLM 配置
    if "llm_config" in config:
        llm = config["llm_config"]
        if "CHAT_MODEL_API" in llm:
            LLMConfig.CHAT_MODEL_API.update(llm["CHAT_MODEL_API"])
        if "MULTIMODAL_CHAT_MODEL_API" in llm:
            LLMConfig.MULTIMODAL_CHAT_MODEL_API.update(llm["MULTIMODAL_CHAT_MODEL_API"])
        LLMConfig.ENABLE_THINKING = llm.get("ENABLE_THINKING", True)
        LLMConfig.ENABLE_SEARCH = llm.get("ENABLE_SEARCH", False)

    return config


def setup_llm(llm_config: LLMConfig = None):
    """
    初始化 LangChain ChatOpenAI 实例。

    Args:
        llm_config: LLM 配置实例。如果为 None，使用全局 LLMConfig。

    Returns:
        ChatOpenAI 实例，或 None 如果配置不完整。
    """
    if llm_config is None:
        llm_config = LLMConfig

    if not llm_config.CHAT_MODEL_API.get("api"):
        print("Warning: LLM API key not configured.")
        return None

    try:
        from langchain_openai import ChatOpenAI
        return ChatOpenAI(
            api_key=llm_config.CHAT_MODEL_API["api"],
            base_url=llm_config.CHAT_MODEL_API["url"],
            model=llm_config.CHAT_MODEL_API["model"],
            temperature=0,
            max_retries=3,
            extra_body={
                "enable_thinking": llm_config.ENABLE_THINKING,
                "enable_search": llm_config.ENABLE_SEARCH,
            },
        )
    except Exception as e:
        print(f"Warning: Failed to initialize LLM: {e}")
        return None

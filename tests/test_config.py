"""Unit tests for configuration system. Run with: pytest tests/test_config.py -v"""

import os
import tempfile
from pathlib import Path

import pytest

from cellhashtag.config.config import (
    CellHashtagConfig,
    LLMConfig,
    load_config,
    PROFILES,
    _deep_merge,
    _set_nested,
)


class TestPydanticModels:
    def test_default_config(self):
        config = CellHashtagConfig()
        assert config.llm.provider == "dashscope"
        assert config.clustering.max_iterations == 3
        assert config.annotation.max_anno_iter == 5
        assert config.lats.search_params.n_iterations == 10

    def test_llm_config_defaults(self):
        llm = LLMConfig()
        assert llm.temperature == 0.3
        assert llm.max_tokens == 2048
        assert llm.enable_thinking is True
        assert llm.enable_search is False

    def test_custom_config(self):
        config = CellHashtagConfig(
            llm=LLMConfig(model="gpt-4o", temperature=0.1),
        )
        assert config.llm.model == "gpt-4o"
        assert config.llm.temperature == 0.1


class TestLoadConfig:
    def test_load_missing_file_uses_defaults(self, tmp_path):
        config = load_config(str(tmp_path / "nonexistent.yaml"))
        assert config.llm.provider == "dashscope"

    def test_load_from_yaml(self, tmp_path):
        yaml_content = """
llm:
  model: "test-model"
  temperature: 0.5
clustering:
  max_iterations: 10
"""
        yaml_file = tmp_path / "test_config.yaml"
        yaml_file.write_text(yaml_content)
        config = load_config(str(yaml_file))
        assert config.llm.model == "test-model"
        assert config.llm.temperature == 0.5
        assert config.clustering.max_iterations == 10

    def test_profile_fast(self, tmp_path):
        config = load_config(str(tmp_path / "none.yaml"), profile="fast")
        assert config.llm.temperature == 0.1
        assert config.clustering.max_iterations == 1
        assert config.annotation.max_anno_iter == 2

    def test_profile_deep(self, tmp_path):
        config = load_config(str(tmp_path / "none.yaml"), profile="deep")
        assert config.llm.temperature == 0.2
        assert config.clustering.max_iterations == 5
        assert config.annotation.max_anno_iter == 8

    def test_runtime_overrides(self, tmp_path):
        config = load_config(
            str(tmp_path / "none.yaml"),
            llm__model="override-model",
            annotation__max_anno_iter=99,
        )
        assert config.llm.model == "override-model"
        assert config.annotation.max_anno_iter == 99


class TestHelperFunctions:
    def test_set_nested(self):
        data = {}
        _set_nested(data, "a.b.c", "value")
        assert data == {"a": {"b": {"c": "value"}}}

    def test_deep_merge(self):
        base = {"a": 1, "b": {"c": 2, "d": 3}}
        override = {"b": {"c": 99}, "e": 5}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": {"c": 99, "d": 3}, "e": 5}

    def test_profiles_exist(self):
        assert "fast" in PROFILES
        assert "default" in PROFILES
        assert "deep" in PROFILES

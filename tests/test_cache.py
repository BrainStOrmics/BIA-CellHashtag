"""Unit tests for EvidenceCache. Run with: pytest tests/test_cache.py -v"""

import time
import tempfile
from pathlib import Path

import pytest

from cellhashtag.core.cache import EvidenceCache


@pytest.fixture
def cache(tmp_path):
    EvidenceCache.CACHE_DIR = tmp_path
    EvidenceCache.DB_PATH = tmp_path / "test_evidence.db"
    return EvidenceCache()


class TestEvidenceCache:
    def test_put_and_get(self, cache):
        cache.put("key1", {"cell_type": "T cell", "confidence": 0.9}, ttl=3600)
        result = cache.get("key1")
        assert result == {"cell_type": "T cell", "confidence": 0.9}

    def test_get_nonexistent_key(self, cache):
        assert cache.get("nonexistent") is None

    def test_ttl_expiration(self, cache):
        cache.put("expire_me", {"data": "old"}, ttl=0.1)
        time.sleep(0.2)
        assert cache.get("expire_me") is None

    def test_clear_expired(self, cache):
        cache.put("keep", {"data": "fresh"}, ttl=3600)
        cache.put("expire", {"data": "stale"}, ttl=0.1)
        time.sleep(0.2)
        cache.clear_expired()
        assert cache.get("keep") is not None
        assert cache.get("expire") is None

    def test_overwrite_key(self, cache):
        cache.put("key", {"version": 1}, ttl=3600)
        cache.put("key", {"version": 2}, ttl=3600)
        result = cache.get("key")
        assert result == {"version": 2}

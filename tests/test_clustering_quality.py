"""Unit tests for clustering quality skills. Run with: pytest tests/test_clustering_quality.py -v"""

import pytest
import scanpy as sc

from cellhashtag.skills.clustering_quality.silhouette_check import silhouette_check
from cellhashtag.skills.clustering_quality.modularity_check import modularity_check


@pytest.fixture
def adata():
    a = sc.read_h5ad("data/example.h5ad")
    sc.tl.leiden(a, resolution=0.5, key_added="leiden")
    return a


class TestSilhouette:
    def test_basic(self, adata):
        result = silhouette_check(adata, "leiden")
        assert "asw" in result
        assert -1 <= result["asw"] <= 1
        assert result["recommendation"] in ("ok", "adjust", "maybe_overcluster")

    def test_per_cluster(self, adata):
        result = silhouette_check(adata, "leiden")
        assert len(result["per_cluster"]) > 0

    def test_missing_embedding(self, adata):
        with pytest.raises(ValueError, match="not found"):
            silhouette_check(adata, "leiden", use_rep="X_nonexistent")


class TestModularity:
    def test_basic(self, adata):
        result = modularity_check(adata, "leiden")
        assert 0 <= result["modularity"] <= 1
        assert result["recommendation"] in ("ok", "adjust")

    def test_missing_connectivities(self, adata):
        del adata.obsp["connectivities"]
        with pytest.raises(ValueError, match="connectivities"):
            modularity_check(adata, "leiden")

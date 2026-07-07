"""Unit tests for format conversion skills. Run with: pytest tests/test_format_conversion.py -v"""

import pytest
import scanpy as sc
import tempfile
from pathlib import Path

from cellhashtag.skills.format_conversion.h5ad_loader import load_h5ad
from cellhashtag.skills.format_conversion.mtx_loader import load_mtx
from cellhashtag.skills.format_conversion.rds_loader import load_rds


class TestH5adLoader:
    def test_load_valid_h5ad(self):
        """Test loading a valid h5ad file."""
        adata = load_h5ad("data/example.h5ad")
        assert adata is not None
        assert adata.n_obs > 0
        assert adata.n_vars > 0

    def test_load_nonexistent_file(self):
        """Test loading a non-existent h5ad file."""
        with pytest.raises(FileNotFoundError):
            load_h5ad("data/nonexistent.h5ad")


class TestMtxLoader:
    def test_load_nonexistent_dir(self):
        """Test loading from non-existent directory."""
        with pytest.raises(FileNotFoundError):
            load_mtx("data/nonexistent_mtx_dir")


class TestRdsLoader:
    def test_load_nonexistent_file(self):
        """Test loading a non-existent RDS file."""
        with pytest.raises(FileNotFoundError):
            load_rds("data/nonexistent.rds")

    def test_unsupported_method(self, tmp_path):
        """Test with unsupported conversion method."""
        # Create a dummy file
        dummy_rds = tmp_path / "test.rds"
        dummy_rds.write_text("dummy")

        with pytest.raises(ImportError, match="rpy2"):
            load_rds(str(dummy_rds), method="anndata2ri")

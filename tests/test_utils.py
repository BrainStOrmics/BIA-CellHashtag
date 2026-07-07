"""Unit tests for I/O utilities. Run with: pytest tests/test_utils.py -v"""

import pytest
import pandas as pd
import scanpy as sc

from cellhashtag.utils.io import (
    save_adata,
    adata_to_pickle,
    adata_from_pickle,
    extract_markers,
    expression_summary,
    infer_omics_type,
    perceive_data,
    df2markdown_table,
)


@pytest.fixture
def adata():
    return sc.read_h5ad("data/example.h5ad")


class TestIO:
    def test_save_and_load(self, adata, tmp_path):
        out = str(tmp_path / "test.h5ad")
        save_adata(adata, out)
        loaded = sc.read_h5ad(out)
        assert loaded.n_obs == adata.n_obs

    def test_pickle_roundtrip(self, adata, tmp_path):
        pkl = str(tmp_path / "test.pkl")
        adata_to_pickle(adata, pkl)
        loaded = adata_from_pickle(pkl)
        assert loaded.n_obs == adata.n_obs


class TestDataProcessing:
    def test_extract_markers(self, adata):
        markers = extract_markers(adata, "leiden", "0", top_n=10)
        assert len(markers) > 0
        assert all(not g.startswith("MT-") for g in markers)

    def test_expression_summary(self, adata):
        summary = expression_summary(adata, ["CD3D", "CD19", "NONEXISTENT"])
        assert "CD3D" in summary or "CD19" in summary
        assert "NONEXISTENT" not in summary

    def test_infer_omics_scrna(self, adata):
        assert infer_omics_type(adata) == "scRNA"

    def test_infer_omics_spatial(self, adata):
        import numpy as np
        adata.obsm["spatial"] = np.zeros((adata.n_obs, 2))
        assert infer_omics_type(adata) == "spRNA"
        del adata.obsm["spatial"]

    def test_perceive_data(self, adata):
        info = perceive_data(adata)
        assert info["n_cells"] == adata.n_obs
        assert info["omics_type"] == "scRNA"


class TestFormatting:
    def test_df2markdown_table(self):
        df = pd.DataFrame({"gene": ["CD3D"], "expr": [1.5]})
        md = df2markdown_table(df)
        assert "| gene | expr |" in md
        assert "| CD3D | 1.5 |" in md

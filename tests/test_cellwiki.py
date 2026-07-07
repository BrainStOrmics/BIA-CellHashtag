"""Unit tests for CellWiki skill. Run with: pytest tests/test_cellwiki.py -v"""

import pytest

from cellhashtag.skills.cellwiki.query_cell_types import (
    query_cell_type,
    list_cell_types,
    query_markers_for_genes,
    _parse_frontmatter,
)


@pytest.fixture
def wiki_dir(tmp_path):
    ct_dir = tmp_path / "cell_types"
    ct_dir.mkdir()
    (ct_dir / "T cell.md").write_text("""---
standard_name: T cell
display_name: T lymphocyte
cl_id: CL:0000084
parent_type: Lymphocyte
---

**Tissues:** blood, lymph node

**Diseases:** lymphoma

| Marker | Type | Evidence | Source |
|--------|------|----------|--------|
| CD3D | surface | protein | CellMarker |
| CD3E | surface | protein | CellMarker |
| CD4 | surface | protein | CellMarker |
""")
    (ct_dir / "B cell.md").write_text("""---
standard_name: B cell
display_name: B lymphocyte
cl_id: CL:0000236
parent_type: Lymphocyte
---

**Tissues:** blood, bone marrow

**Diseases:** myeloma

| Marker | Type | Evidence | Source |
|--------|------|----------|--------|
| CD19 | surface | protein | CellMarker |
| MS4A1 | surface | protein | CellMarker |
""")
    return tmp_path


class TestCellWiki:
    def test_list_cell_types(self, wiki_dir):
        types = list_cell_types(wiki_dir)
        assert "T cell" in types
        assert "B cell" in types

    def test_query_cell_type(self, wiki_dir):
        info = query_cell_type("T cell", wiki_dir)
        assert info is not None
        assert info["cl_id"] == "CL:0000084"
        assert len(info["markers"]) == 3
        assert "blood" in info["tissues"]

    def test_query_nonexistent(self, wiki_dir):
        assert query_cell_type("Nonexistent", wiki_dir) is None

    def test_query_markers_for_genes(self, wiki_dir):
        results = query_markers_for_genes(["CD3D", "CD4", "ACTB"], wiki_dir)
        assert len(results) > 0
        assert results[0]["cell_type"] == "T lymphocyte"
        assert "CD3D" in results[0]["overlap"]

    def test_parse_frontmatter(self):
        content = "---\nname: Test\nvalue: 42\n---\nBody"
        fm = _parse_frontmatter(content)
        assert fm["name"] == "Test"

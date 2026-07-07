---
name: cellwiki
description: Query CellWiki knowledge base for marker genes, cell type ontologies, and tissue-specific cell type constraints.
allowed-tools: Read, Bash, Glob
compatibility: Cell annotation evidence gathering
---

# CellWiki Knowledge Base Query

## When to Use

During cell type annotation to gather evidence for/against candidate hypotheses. Query CellWiki for known markers, cell type definitions, and biological constraints.

## Data Location

CellWiki wiki files are at `~/CellWiki/wiki/` with:
- `cell_types/*.md` — One file per cell type with frontmatter and markers table
- `marker_genes/` — (future) gene-centric marker database
- `tissues/` — Tissue-specific cell type compositions
- `contradictions.md` — Known biologically incompatible marker combinations

## Query Functions

### `query_cell_type(name, wiki_dir) -> dict`

Reads `~/CellWiki/wiki/cell_types/{name}.md` and returns:
```python
{
    "name": str,           # standard_name from frontmatter
    "display_name": str,
    "cl_id": str,          # Cell Ontology ID
    "parent_type": str,
    "markers": [{"gene": str, "type": str, "evidence": str}],
    "tissues": list[str],
    "diseases": list[str],
}
```

### `list_cell_types(wiki_dir) -> list[str]`

Lists all available cell type files in `cell_types/`.

### `query_markers_for_genes(genes, wiki_dir) -> list[dict]`

Scans all cell type files, returns those whose markers overlap with the given gene list. Each result:
```python
{"cell_type": str, "overlap": list[str], "jaccard": float, "markers_all": list[str]}
```

### `get_tissue_cell_types(tissue, wiki_dir) -> list[str]`

Reads `~/CellWiki/wiki/tissues/` to find cell types plausible for a given tissue.

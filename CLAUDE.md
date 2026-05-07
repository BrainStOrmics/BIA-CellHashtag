# BIA-CellHashtag (Cell#) — LLM-Based Single-Cell Annotation Agent

## Project Overview

Automated cell type annotation tool for scRNA-seq and spatial RNA-seq data using LLMs. Uses self-criticism mechanisms and web scraping to achieve expert-level annotation quality. Built on LangGraph.

## Architecture

```
src/cellhashtag/
  _Agent.py         # Main agent class
  _prompts.py       # Prompt templates
  graph/
    cellhashtag.py  # Main annotation graph
    annotator.py    # Annotation logic
  config/
    config.py       # Configuration
    config.yaml     # User config
  utils/
    io.py           # Data I/O
    utilities.py    # Helper functions
    setup.py        # Setup utilities
  data/
    example.h5ad    # Example AnnData
    example_cellmarkers.csv  # Marker reference
```

## Dependencies

- Python (check pyproject.toml if present)
- LangGraph for workflow orchestration
- scanpy / AnnData for single-cell data handling
- LLM API (OpenAI-compatible or similar)
- Web scraping for external knowledge (CellMarker, PanglaoDB)

## LangGraph Workflow

1. **Data Loading**: Load AnnData (.h5ad) with count matrix
2. **QC & Preprocessing**: Filter cells, normalize, cluster
3. **Marker Identification**: Find cluster-specific marker genes
4. **LLM Annotation**: LLM annotates cell types based on markers
5. **Self-Criticism**: LLM evaluates its own annotations
6. **Revision**: If critique fails, revise with additional context
7. **Output**: Final annotations with confidence scores

## Data Format

Input: AnnData (.h5ad) with:
- `.X`: Count matrix (cells × genes)
- `.obs`: Cell metadata
- `.var`: Gene metadata
- `.obsm`: Embeddings (UMAP, PCA)

Output: Cell type annotations added to `.obs['cell_type']`

## Code Standards

- Type hints on all public functions
- LangGraph state only contains lightweight data (paths, summaries)
- Never put full AnnData objects in graph state
- Use existing prompt templates in `_prompts.py`
- Configuration via `config/` — don't hardcode model endpoints
- Self-criticism loops must have max iteration limits

## Common Patterns

- Marker gene extraction → LLM prompt → annotation → critique loop
- Web scraping for external marker databases (CellMarker, PanglaoDB)
- Confidence scoring for all annotations
- Support for both scRNA-seq and spatial RNA-seq
- Flexible data compatibility (scanpy, AnnData formats)

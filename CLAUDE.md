# BIA-CellHashtag (Cell#) — LLM-Based Single-Cell Annotation Agent

## Project Overview

Automated cell type annotation tool for scRNA-seq and spatial RNA-seq data using LLMs. Uses self-criticism mechanisms and web scraping to achieve expert-level annotation quality. Built on LangGraph.

## Architecture

```
src/cellhashtag/
  _Agent.py              # Main agent class
  _prompts.py            # Prompt templates
  graph/
    orchestrator.py      # Top-level workflow: clustering → annotation → summary
    high_hierarchy/
      __init__.py        # Exports create_graph()
      graph.py           # High hierarchy state machine
      nodes.py           # DataRep, Estimate, Knowledge, Reflection
    low_hierarchy/
      graph.py
      nodes.py           # DataRep, Harness
    tree_search/
      graph.py
      nodes.py           # Search, Memory
  config/
    config.py            # Configuration
    config.yaml          # User config
  utils/
    io.py                # Data I/O
    utilities.py         # Helper functions
    setup.py             # Setup utilities
  data/
    example.h5ad         # Example AnnData
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
2. **Adaptive Clustering**: Filter cells, normalize, auto-determine optimal clustering resolution
3. **High Hierarchy Annotation**: DataRep → Estimate → Knowledge → Reflection (推理循环)
4. **Low Hierarchy Swarm** (if Reflection triggers): DataRep → Harness → Search  Memory (并行树搜索)
5. **Results Summary**: Final annotations with confidence scores and report

## Data Format

Input: AnnData (.h5ad) with:
- `.X`: Count matrix (cells × genes)
- `.obs`: Cell metadata
- `.var`: Gene metadata
- `.obsm`: Embeddings (UMAP, PCA)

Output: Cell type annotations added to `.obs['cell_type']`

## Code Standards

- Type hints on all public functions
- LangGraph state only contains lightweight data (paths, summaries) — never full AnnData objects
- Each subgraph lives in its own directory (`graph/<name>/`) with `graph.py` + `nodes.py`
- Node functions are scoped to their parent graph only — no cross-graph sharing
- Subgraphs nest unidirectionally: `orchestrator → high_hierarchy → tree_search`
- Prompts via `prompts/` directory, not hardcoded in Python
- Configuration via `config/` — don't hardcode model endpoints
- All annotation/criticism loops must have `max_anno_iter` limits

## Common Patterns

- High hierarchy: marker gene extraction → LLM hypothesis → external knowledge → reflection decision
- Low hierarchy swarm: constraint-driven tree search with memory for sub-cluster exploration
- Web scraping for external marker databases (CellMarker, PanglaoDB)
- Confidence scoring for all annotations
- Support for both scRNA-seq and spatial RNA-seq

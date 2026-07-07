# BIA-CellHashtag (Cell#) v3.1 — LLM-Based Single-Cell Annotation Agent

Automated cell type annotation for scRNA-seq and spatial RNA-seq data using LLMs. Built on LangGraph + DeepAgents with LATS tree search, iterative clustering quality assessment, and multi-omics support.

> **v3.1 adds 8 new skills (25 exports): QC, DEG, plotting, batch correction, spatial loaders, trajectory inference.**

## Features

- **LATS Tree Search**: MCTS-based hypothesis exploration with evidence caching (SQLite + TTL)
- **DeepAgents Integration**: Multi-agent annotation with self-criticism and knowledge retrieval
- **Iterative Clustering**: Leiden with silhouette + modularity quality checks and automatic resolution tuning
- **Multi-format Input**: `.h5ad`, `.rds` (Seurat/SCE), 10x `.mtx`, Visium, MERFISH, Stereo-seq
- **Multi-omics Support**: scRNA-seq, snRNA-seq, spatial RNA-seq, ATAC-seq, CITE-seq
- **CellWiki Integration**: Structured knowledge base for marker genes and cell types
- **Pydantic Configuration**: Type-safe config with YAML + env var + runtime override layers
- **Rich Output**: Annotated AnnData, CSV table, Markdown report, UMAP visualization, HITL review
- **Quality Control**: Mitochondrial/ribosomal metrics, cell filtering, doublet detection
- **Batch Correction**: Harmony, BBKNN, ComBat integration
- **Trajectory Inference**: PAGA graph-based developmental trajectories

## Installation

```bash
pip install -e .
```

### Dependencies

- Python 3.11+
- `langgraph`, `langchain-openai` — workflow orchestration and LLM
- `scanpy`, `anndata` — single-cell data handling
- `scikit-learn` — silhouette score
- `python-igraph` — modularity computation
- `pyyaml` — configuration loading
- `pandas`, `numpy` — data manipulation

## Quick Start

```python
from cellhashtag import CellHashtagAgent

# Create agent (loads config from src/cellhashtag/config/config.yaml)
agent = CellHashtagAgent(profile="fast")  # or "default", "deep"

# Run annotation pipeline
result = agent.run(
    input_path="data/example.h5ad",
    output_dir="output",
    cluster_key="leiden",
)

print(f"Status: {result['status']}")
print(f"Clusters: {result['n_clusters']}")
```

### Configuration

Edit `src/cellhashtag/config/config.yaml` or override at runtime:

```python
agent = CellHashtagAgent(
    profile="default",
    llm__model="qwen-plus",
    annotation__max_anno_iter=3,
)
```

See `.research_harness_system/QUICKSTART.md` for full setup guide.

## Configuration

Copy and edit the config file:

```bash
cp src/cellhashtag/config/config.yaml.example src/cellhashtag/config/config.yaml
```

Edit `config.yaml` to set your LLM API key, URL, and model:

```yaml
llm_config:
  CHAT_MODEL_API:
    api: "your-api-key"
    url: "https://your-api-endpoint"
    model: "your-model-name"
    type: openai
  ENABLE_THINKING: true
  ENABLE_SEARCH: false
```

## Architecture

v3.0 uses a flat 4-node LangGraph orchestrator:

```
perception → clustering → annotation → output
```

- **Perception**: Load h5ad, infer omics type, extract top markers
- **Clustering**: Iterative Leiden with silhouette + modularity quality checks, automatic resolution tuning
- **Annotation**: DeepAgents + LATS MCTS tree search per cluster (with error handling)
- **Output**: Save annotated h5ad + CSV table + Markdown report + UMAP plots + HITL review

### Project Structure

```
src/cellhashtag/
├── agent.py              # CellHashtagAgent entry point
├── config/
│   ├── config.py         # Pydantic models + load_config()
│   └── config.yaml       # User configuration
├── core/
│   ├── mcts.py           # Pure MCTS algorithm (no LLM/IO)
│   ├── lats_loop.py      # LATS outer loop: wires MCTS <-> DeepAgents
│   └── cache.py          # SQLite evidence cache with TTL
├── graphs/
│   ├── orchestrator.py   # 4-node flat graph
│   └── clustering.py     # Iterative clustering subgraph
├── skills/
│   ├── cellwiki/            # CellWiki knowledge queries
│   ├── clustering_quality/  # Silhouette, modularity, phiclust
│   ├── clustering_degs/     # find_markers, consensus, resolution scan
│   ├── format_conversion/   # h5ad, rds, mtx, SCE loaders
│   ├── quality_control/     # QC metrics, filtering, doublets
│   ├── batch_correction/    # Harmony, BBKNN, ComBat
│   ├── trajectory/          # PAGA trajectory inference
│   ├── spatial/             # Visium, MERFISH, Stereo-seq, generic
│   ├── lats_search/         # MCTS tests
│   └── plotting/
│       ├── embeddings/      # UMAP/t-SNE (Okabe-Ito, PDF+PNG)
│       ├── gene_expression/ # violin, dot plots
│       └── heatmaps/        # marker heatmaps with dendrogram
├── utils/
│   ├── io.py             # AnnData I/O, marker extraction
│   └── hitl.py           # Human-in-the-Loop review
└── prompts/
    ├── annotation.md                # Annotation system prompt
    ├── annotation_hypothesis.md     # Hypothesis generation prompt
    ├── annotation_evaluation.md     # Adversarial evaluation prompt
    ├── annotation_task.md           # Per-cluster task template
    └── lats_search.md               # LATS subagent prompt
```

See `.research_harness_system/ARCHITECTURE.md` for full details.

## Development

```bash
# Install dependencies (using BIA-dev conda environment)
conda activate BIA-dev

# Run MCTS unit tests
python -m pytest src/cellhashtag/skills/lats_search/test_mcts_core.py -v

# Interactive component tests (Jupyter notebook, 14 incremental cells)
jupyter notebook local_tests/test_lats_components.ipynb
```

## Documentation

- **Quick Start**: `.research_harness_system/QUICKSTART.md`
- **Architecture**: `.research_harness_system/ARCHITECTURE.md`
- **Demo Notebook**: `examples/DEMO.ipynb`
- **Usage Example**: `examples/basic_usage.py`

## License

[MIT License](LICENSE)

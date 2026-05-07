# BIA-CellHashtag (Cell#) v2.0 — LLM-Based Single-Cell Annotation Agent

Automated cell type annotation for scRNA-seq and spatial RNA-seq data using LLMs. Built on LangGraph with built-in clustering quality assessment, multi-omics support, and a self-criticism annotation loop.

## Features

- **LLM Self-Criticism**: Annotation → Critique → Revision loop for expert-level quality
- **Built-in Clustering**: Leiden clustering with multi-metric quality assessment (Silhouette, Modularity, Phiclust) and automatic resolution tuning
- **Multi-format Input**: `.h5ad`, `.rds` (Seurat), 10x `.mtx`
- **Multi-omics Support**: scRNA-seq, snRNA-seq, spatial RNA-seq, ATAC-seq, CITE-seq, Metabolomics
- **Dual-layer Annotation**: Layer 1 (Cell Type) + Layer 2 (Subpopulation / State / Function / Spatial)
- **CellWiki Integration**: Structured knowledge base for marker genes and cell types
- **HITL (Human-in-the-Loop)**: Interactive checkpoints at clustering quality分歧, metadata verification, and annotation review
- **Rich Output**: Updated AnnData, CSV annotation table, Markdown report with UMAP visualization

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
from cellhashtag import CellHashtagAgent, setup_llm

# Initialize LLM (configure in config/config.yaml or pass directly)
llm = setup_llm()

# Create agent
agent = CellHashtagAgent()

# Run annotation
result = agent.run(
    input_path="data/pbmc3k_raw.h5ad",
    llm=llm,
    input_format="h5ad",
    cluster_key="leiden",
    omics_type="scRNA",
    cell_marker_df_dir="data/pbmc3k_markers.csv",
    max_iterations=3,      # clustering iterations
    max_anno_iter=5,       # annotation self-criticism iterations
    output_dir="output",
)

print(result["cell_type_annotations"])
```

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

### Pipeline

```
Data Perception → Metadata Extraction → Clustering → Quality Assessment →
Feature Extraction → Knowledge Match → Layer 1 Annotation → Layer 2 Annotation →
Confidence Assessment → Output Generation
```

### Clustering Quality Loop

The clustering node runs an iterative loop:

1. Run Leiden clustering at current resolution
2. Assess quality with Silhouette, Modularity, and Phiclust
3. Decide: proceed / adjust resolution / trigger HITL
4. Repeat up to `max_iterations`

### Annotation Self-Criticism Loop

For each cluster:

1. LLM annotates cell type based on marker genes and metadata
2. Critic LLM evaluates the annotation against 5 criteria
3. If disapproved, revise with critique feedback
4. Repeat up to `max_anno_iter`

### Project Structure

```
src/cellhashtag/
├── agent.py                          # Main agent class
├── state.py                          # LangGraph state definition
├── config/
│   ├── config.py                     # Configuration loader
│   └── config.yaml                   # User config
├── graphs/
│   ├── orchestrator.py               # Main pipeline orchestrator
│   └── clustering.py                 # Clustering subgraph
├── nodes/
│   ├── data_perception.py            # Input data perception
│   ├── metadata_node.py              # Metadata extraction
│   └── clustering_node.py            # Clustering + quality decision
├── skills/
│   ├── clustering_quality/           # Quality assessment tools
│   │   ├── silhouette_check.py       # Silhouette score
│   │   ├── modularity_check.py       # Modularity score
│   │   └── phiclust_check.py         # Phiclust (optional)
│   └── format_conversion/            # Multi-format loaders
│       ├── h5ad_loader.py
│       ├── rds_loader.py
│       └── mtx_loader.py
├── prompts/
│   ├── annotation_layer1.md          # Layer 1 cell type prompt
│   ├── critic.md                     # Self-criticism prompt
│   └── functional.md                 # Layer 2 functional prompt
├── utils/
│   ├── cellwiki_client.py            # CellWiki knowledge interface
│   ├── hitl.py                       # Human-in-the-loop CLI
│   └── utilities.py                  # Helper functions
└── data/
    ├── example.h5ad                  # Example AnnData
    └── example_cellmarkers.csv       # Marker reference
```

## Input / Output

### Input

AnnData (`.h5ad`) with:
- `.X`: Count matrix (cells x genes)
- `.obs`: Cell metadata
- `.var`: Gene names
- `.obsm`: Embeddings (UMAP, PCA) — optional, computed if missing

### Output

- **Annotated AnnData** (`.pkl`): `adata.obs["Cell#"]` contains cell type annotations
- **CSV Table**: Cluster-to-cell-type mapping with confidence scores
- **Markdown Report**: Annotation summary with reasoning

## Compared to v1 (MVP)

| Feature | v1 (MVP) | v2.0 |
|---|---|---|
| Clustering | External, fixed | Built-in, quality-assessed, iterative |
| Multi-omics | scRNA-seq only | scRNA / snRNA / spRNA / ATAC / CITE / Metab |
| Knowledge | LLM internal | CellWiki integration |
| Confidence | None | Multi-dimensional scoring |
| HITL | No | Clustering + Metadata + Annotation review |
| Output | UMAP only | AnnData + CSV + Markdown report |
| State design | AnnData in state | File paths only, LangGraph-safe |

## License

[MIT License](LICENSE)

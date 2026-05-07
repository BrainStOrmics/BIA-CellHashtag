# BIA-CellHashtag v2.0 重构方案

> 从 MVP 到生产级多组学细胞注释 Agent
> 基于 CellWiki 知识底座 + LangGraph 状态机 + 工具 Skill 集成
> 生成日期: 2026-05-05

---

## 一、现有 MVP 分析

### 1.1 架构概览

```
MainGraph:
  START → get_metadata → get_targetcell → [scrape_markers?] → get_expression → annotation_subgraph → normalize_celltype → END

SubGraph (per cluster):
  START → annotation → critic → [finalize | reflex→annotation] → [functional | END]
```

### 1.2 优点

- **LLM 自批评机制**：annotation → critic → reflex 循环提高注释质量
- **Map-Reduce 模式**：每个 cluster 独立注释，最后 normalize 统一
- **Web 搜索增强**：通过 Tavily 获取外部 marker 信息
- **LangGraph 状态机**：使用 checkpointer 支持断点续跑

### 1.3 待改进

| 问题 | 现状 | v2.0 方案 |
|------|------|-----------|
| 聚类质量评估 | ❌ 无，聚类在 Agent 外进行 | ✅ Clustering Quality Assessment Skill（Phiclust + Silhouette + Clustree） |
| 聚类迭代 | ❌ 无，固定 resolution | ✅ Clustering ↔ Feature Extraction 迭代循环，自动调整 resolution |
| 多组学支持 | ❌ 仅 scRNA-seq | ✅ 多组学路由子图（scRNA / snRNA / spRNA / ATAC / CITE-seq / Metabolism） |
| 知识底座 | ❌ 无，依赖 LLM 内部知识 | ✅ 集成 CellWiki 作为结构化知识查询引擎 |
| Metadata 来源 | ⚠️ 仅 adata.obs | ✅ obs + 用户输入表格 + HITL 核对 |
| 双层注释 | ⚠️ 有 functional 但无 subpopulation 层次 | ✅ Layer 1: Cell Type + Layer 2: 亚群/状态/功能/空间 |
| HITL 交互 | ❌ 无 | ✅ 聚类质量分歧 + Metadata 核对 + 注释审查 |
| 数据格式 | ⚠️ 仅 AnnData | ✅ 多格式输入 → 统一 AnnData 内部处理 |
| 输出 | ⚠️ 仅 UMAP 图 | ✅ 更新 AnnData.obs + 注释表格 + Markdown 报告（UMAP + heatmap） |
| 置信度评估 | ⚠️ 无结构化 | ✅ 多维度置信度（marker_match / knowledge_consistency / cluster_quality / multi_omics_agreement） |
| 证据层级 | ❌ 无 | ✅ 继承 CellWiki 的 1-5 证据层级 |

---

## 二、v2.0 架构设计

### 2.1 工作流概览

```
┌─────────────────────────────────────────────────────────────────┐
│                    CellHashtag v2.0 Orchestrator                  │
│                                                                   │
│  1. Data Perception (格式转换 + QC 验证)                         │
│         │                                                         │
│         ▼                                                         │
│  2. Metadata Extraction (obs + 用户输入 + HITL 核对)              │
│         │                                                         │
│         ▼                                                         │
│  3. Clustering ──────────────────────┐                           │
│         │                              │                           │
│         ▼                              │                           │
│  4. Clustering Quality Assessment ◄───┘ 迭代调整 resolution       │
│     (Phiclust + Silhouette + Clustree)   (max N 次)               │
│         │                                                         │
│         ▼                                                         │
│  5. Feature Extraction (按组学路由)                               │
│     scRNA ─→ DE analysis                                          │
│     ATAC  ─→ peak-to-gene linking                                 │
│     CITE  ─→ protein DE                                           │
│     spRNA ─→ spatial-aware DE                                     │
│         │                                                         │
│         ▼                                                         │
│  6. Knowledge Match (查询 CellWiki)                               │
│     marker_genes + cell_types + trajectories + metabolites        │
│         │                                                         │
│         ▼                                                         │
│  7. Layer 1 Annotation (Cell Type)                                │
│         │                                                         │
│         ▼                                                         │
│  8. Layer 2 Annotation (并行子图)                                 │
│     ├─ Subpopulation (亚群细分)                                    │
│     ├─ State (状态: activated / proliferating / ...)              │
│     ├─ Function (功能: cytotoxic / immunosuppressive / ...)       │
│     └─ Spatial (空间: tumor_core / invasive_margin / ...)         │
│         │                                                         │
│         ▼                                                         │
│  9. Confidence Assessment (多维度)                                │
│         │                                                         │
│         ▼                                                         │
│  10. Output Generation                                            │
│     - Updated AnnData (.obs 列)                                   │
│     - 注释表格 (cell_id → hashtags)                               │
│     - Markdown 报告 (UMAP + heatmap + 解释)                       │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 LangGraph 状态机设计

#### CellHashtagState

```python
class CellHashtagState(TypedDict):
    # ===== 输入 =====
    input_path: str                      # 输入文件路径
    input_format: str                    # h5ad / rds / mtx / csv
    metadata_path: str                   # 可选：额外 metadata CSV
    omics_type: str                      # scRNA / snRNA / spRNA / ATAC / CITE / Metab / multi
    
    # ===== Data Perception =====
    perceived_info: dict                 # 物种、平台、细胞数、基因数等
    adata_path: str                      # 标准化后的 AnnData 路径（不在 state 中存对象）
    
    # ===== Metadata =====
    metadata_dict: dict                  # {cluster_id: {metadata_key: value}}
    metadata_summary: str               # Metadata 摘要
    metadata_confirmed: bool             # HITL 是否已确认
    
    # ===== Clustering =====
    cluster_key: str                     # adata.obs 中的聚类列名
    clustering_params: dict              # {resolution, n_neighbors, n_pcs, ...}
    clustering_result: dict              # {cluster_labels, n_clusters}
    clustering_history: list[dict]       # [{iteration, params, quality, decision}]
    clustering_quality: dict             # 质量评估结果
    iteration_count: int                 # 当前迭代次数
    max_iterations: int                  # 最大迭代次数
    
    # ===== Feature Extraction =====
    features: dict                       # {omics_type: {cluster_id: marker_list}}
    
    # ===== Knowledge Match =====
    matched_types: dict                  # {cluster_id: [(cell_type, score), ...]}
    contradictions: list[dict]           # 标记矛盾列表
    evidence_tiers: dict                 # {cluster_id: evidence_tier}
    
    # ===== Layer 1: Cell Type =====
    cell_type_annotations: dict          # {cluster_id: cell_type}
    layer1_confidence: dict              # {cluster_id: confidence_score}
    
    # ===== Layer 2: Sub-annotation =====
    subtype_annotations: dict            # {cluster_id: [subtypes]}
    state_annotations: dict              # {cluster_id: [states]}
    function_annotations: dict           # {cluster_id: [functions]}
    spatial_annotations: dict            # {cluster_id: [spatial_zones]}
    layer2_confidence: dict              # {cluster_id: {subtype: score, ...}}
    
    # ===== Confidence =====
    confidence: dict                     # 多维度置信度
    
    # ===== Output =====
    hashtags: dict                       # {cluster_id: ["#CD4_T_cell", "#regulatory", ...]}
    report_path: str                     # Markdown 报告路径
    annotation_table_path: str           # CSV 表格路径
    output_adata_path: str               # 输出 AnnData 路径
    
    # ===== Status =====
    status: str                          # perception → metadata → clustering → ... → done
    errors: list[str]                    # 错误列表
    interrupt_message: str               # HITL 消息
```

### 2.3 工具 Skill 集成

#### Clustering Quality Assessment Skill

| 工具 | 功能 | 决策权重 |
|------|------|---------|
| **Phiclust** | 基于随机矩阵理论，检测 cluster 内隐藏亚群 | 0.4 |
| **Silhouette** | 轮廓系数，衡量 cluster 内聚 vs 间分离 | 0.3 |
| **Modularity** | Leiden/Louvain 模块度 | 0.2 |
| **Clustree** | 多 resolution 聚类树可视化 | 0.1 |

**决策逻辑**：
- 多数工具认为需要调整 → 自动调整 resolution 并重试
- 工具意见分歧 → 触发 HITL
- 多数工具认为 OK → 继续下一步

#### 数据格式转换 Skill

| 输入格式 | 转换工具 | 说明 |
|---------|---------|------|
| .h5ad | 直接加载 | AnnData 原生格式 |
| .rds (Seurat) | `anndata2ri` / `zellkonverter` | R → Python 转换 |
| 10x mtx | `scanpy.read_10x_mtx()` | 原生支持 |
| .loom | `scanpy.read_loom()` | 原生支持 |
| CSV/TSV | `scanpy.read_text()` | 需指定分隔符 |

### 2.4 HITL 交互设计

**交互时机**：
1. **Metadata 核对**：Agent 感知 Metadata 后，与用户确认是否正确
2. **聚类质量分歧**：工具评估意见不一致时，让用户拍板
3. **注释审查**（可选）：关键 cluster 注释完成后，允许用户手动调整

**交互方式**：
- CLI 交互式提示
- Markdown 报告 + 可视化
- 用户通过配置文件覆盖自动决策

---

## 三、目录结构

```
src/cellhashtag/
├── __init__.py
├── agent.py                    # 主 Agent 类（原 _Agent.py 重写）
├── state.py                    # LangGraph 状态定义
├── config/
│   ├── __init__.py
│   ├── config.py               # 配置加载
│   └── config.yaml             # 用户配置
│
├── graphs/                     # LangGraph 子图
│   ├── __init__.py
│   ├── orchestrator.py         # 主图编排
│   ├── clustering.py           # 聚类 + 质量评估子图
│   ├── feature_extraction.py   # 特征提取（按组学路由）
│   ├── knowledge_match.py      # CellWiki 查询匹配
│   ├── annotation_layer1.py    # Layer 1: Cell Type
│   ├── annotation_layer2.py    # Layer 2: 亚群/状态/功能/空间
│   ├── confidence.py           # 置信度评估
│   └── output.py               # 输出生成
│
├── skills/                     # 工具 Skill
│   ├── __init__.py
│   ├── clustering_quality/     # 聚类质量评估
│   │   ├── __init__.py
│   │   ├── phiclust_check.py   # Phiclust
│   │   ├── silhouette_check.py # Silhouette
│   │   ├── modularity_check.py # Modularity
│   │   └── clustree_check.py   # Clustree
│   ├── format_conversion/      # 数据格式转换
│   │   ├── __init__.py
│   │   ├── h5ad_loader.py
│   │   ├── rds_loader.py
│   │   └── mtx_loader.py
│   └── metadata_extraction/    # Metadata 提取
│       ├── __init__.py
│       └── obs_parser.py
│
├── nodes/                      # LangGraph 节点实现
│   ├── __init__.py
│   ├── data_perception.py
│   ├── metadata_node.py
│   ├── clustering_node.py
│   ├── feature_extraction_node.py
│   ├── knowledge_match_node.py
│   ├── annotation_nodes.py
│   ├── confidence_node.py
│   └── output_node.py
│
├── prompts/                    # Prompt 模板
│   ├── __init__.py
│   ├── annotation_layer1.md
│   ├── annotation_layer2.md
│   ├── critic.md
│   └── functional.md
│
├── utils/                      # 工具函数
│   ├── __init__.py
│   ├── hitl.py                 # HITL 交互
│   ├── plotting.py             # 可视化（UMAP + heatmap）
│   ├── report.py               # Markdown 报告生成
│   └── cellwiki_client.py      # CellWiki 查询接口
│
└── data/                       # 示例数据
    ├── example.h5ad
    └── example_cellmarkers.csv
```

---

## 四、与现有 MVP 的兼容性

### 4.1 保留的设计

- ✅ LangGraph 状态机架构
- ✅ LLM 自批评机制（annotation → critic → reflex）
- ✅ Map-Reduce 模式（每个 cluster 独立处理）
- ✅ Web 搜索增强（可选）
- ✅ 配置文件系统（config.yaml）

### 4.2 新增的设计

- ✅ Clustering 内置 + 质量评估迭代
- ✅ 多组学路由子图
- ✅ CellWiki 集成
- ✅ 双层注释体系
- ✅ HITL 交互
- ✅ 多维度置信度
- ✅ 丰富的输出（AnnData + CSV + Markdown 报告）

### 4.3 废弃的设计

- ❌ 外部聚类（现在内置）
- ❌ 简单 UMAP 输出（现在多维度输出）
- ❌ 无结构化置信度
- ❌ pickle 序列化 AnnData（现在用文件路径）

---

## 五、实施路线图

### Phase 0: 基础设施
- [ ] 创建新目录结构
- [ ] 定义 CellHashtagState
- [ ] 实现配置加载

### Phase 1: 核心工作流
- [ ] Data Perception 节点
- [ ] Metadata 提取 + HITL
- [ ] Clustering 节点
- [ ] Clustering Quality Assessment Skill

### Phase 2: 知识集成
- [ ] CellWiki 查询接口
- [ ] Knowledge Match 节点
- [ ] 格式转换 Skill

### Phase 3: 注释引擎
- [ ] Layer 1 Annotation（Cell Type）
- [ ] Layer 2 Annotation（亚群/状态/功能/空间）
- [ ] LLM 自批评机制

### Phase 4: 输出与报告
- [ ] 置信度评估
- [ ] 注释表格生成
- [ ] Markdown 报告 + 可视化

### Phase 5: 多组学扩展
- [ ] 组学路由子图
- [ ] ATAC / CITE-seq / Metabolism 支持

---

*本文档为 BIA-CellHashtag v2.0 重构方案，基于现有 MVP 进行架构升级。*
*核心原则：保留 LLM 自批评和 Map-Reduce 的优点，引入聚类质量评估、多组学支持和 CellWiki 集成。*

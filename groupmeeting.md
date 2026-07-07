# 🧬 CellHashtag (Cell#) Project Report — Group Meeting

## 1. Project Overview

Cell# 是一个基于 LLM 的单细胞自动注释工具，利用大语言模型和 MCTS 树搜索实现专家级细胞类型注释质量。输入单细胞测序数据（AnnData），输出每个 cluster 的细胞类型标注及其置信度。

> 🎯 **核心定位**: 将人工注释从数天缩短到数十分钟，同时保持专家级注释质量。

---

## 2. 版本演进历史

### 🔨 v0.x → v1.x (Initial Development)

- 基础 pipeline 实现：数据加载 → 聚类 → LLM 注释
- 手写 annotation-critic 循环
- 硬编码配置和 prompt

### 🚀 v2.0 (Major Release)

- 双层演进状态机架构：High Hierarchy（推理循环）+ Low Hierarchy（并行搜索）
- 引入 Web scraping（CellMarker、PanglaoDB）获取外部知识
- 自我批评机制：annotation → critic → reflection 迭代
- 配置系统初版：`config.yaml` 手动解析

### 🔄 v3.0 (Architecture Rewrite) — 2026-05-26

**Breaking changes from v2.0**: 删除 high_hierarchy、low_hierarchy、nodes/、state.py

- **Flat 4-node LangGraph orchestrator**: perception → clustering → annotation → output
- **DeepAgents 替代手写 annotation-critic 循环**，大幅降低 orchestrator 复杂度
- **LATS MCTS 作为 DeepAgents subagent** 运行，保留纯算法模块（core/mcts.py）
- **Pydantic 配置模型**替代 dict，支持 profile（fast/default/deep）和环境变量覆盖
- **SQLite EvidenceCache** 替代 pickle 缓存，支持 TTL 过期和并发安全
- **CellWiki 本地知识库**替代在线 web scraping，提供可靠的标记基因证据

### 🌟 v3.1 (Skills Expansion) — Current

- **8 个新 skills**（25 个函数导出），覆盖 QC、DEG、绘图、批次校正、空间组学、轨迹推断
- **Prompt 外部化**：3 个 prompt 从代码中提取到 `prompts/*.md`
- **Markdown 报告 + UMAP 可视化**输出
- **HITL（Human-in-the-Loop）重新实现**：低置信度注释交互式审核
- **注解错误处理**：单个 cluster 失败不影响全局流程
- **全量技术债清零** ✅

---

## 3. 当前架构 (v3.1)

### 3.1 核心 Workflow — 4 Node 流水线

```
┌─────────────┐     ┌──────────────┐     ┌───────────────┐     ┌──────────┐
│ perception  │────▶│ clustering   │────▶│ annotation    │────▶│ output   │
│ 数据感知     │     │ 聚类优化      │     │ LLM 注释       │     │ 结果输出  │
└─────────────┘     └──────────────┘     └───────────────┘     └──────────┘
```

| Node | 输入 | 输出 | 关键操作 |
|------|------|------|----------|
| **perception** 📊 | `.h5ad` 文件路径 | 聚类 key、标记基因字典、omics type | 加载 AnnData → 过滤/标准化/PCA → Leiden 聚类 → 提取每个 cluster top 标记基因 |
| **clustering** 🔗 | 初始聚类结果 | 最优聚类 + 质量评估 | 多工具投票（silhouette + modularity + phiclust）→ 自动调节 resolution → 迭代优化 |
| **annotation** 🤖 | 聚类结果 + 标记基因 | 每个 cluster 的细胞类型标注 | DeepAgents 主 agent + LATS MCTS subagent + CellWiki 知识检索 |
| **output** 📦 | 标注结果 | annotated.h5ad + CSV + Markdown report + UMAP | 写入 AnnData（`Cell#`, `Cell#_confidence` 列）+ HITL 审核 + 可视化 |

### 3.2 Annotation 阶段详解

```
DeepAgents (主 agent) 🧠
  ├── System Prompt: 细胞类型注释专家角色
  ├── Skills:
  │   ├── cellwiki 📚 (标记基因查询、细胞类型本体)
  │   └── clustering_quality 🔍 (质量评估)
  └── SubAgent: lats-search 🔎
        ├── MCTS 纯算法模块 (core/mcts.py)
        │   └── UCB1 selection → expand → validate → evaluate → backpropagate
        └── CellWiki 知识检索 (Jaccard 相似度匹配)
```

每个 cluster 串行调用 `deep_agent.invoke()`，带 try/except 错误处理：
- ✅ 成功 → 返回 JSON `{cell_type, confidence, reasoning}`
- ❌ 失败 → 标记为 "Unknown"，记录错误信息，继续处理后续 cluster

---

## 4. Skills 系统 (25 个函数)

### 4.1 完整 Skills 列表

| 分类 | Skill | 函数 | 说明 |
|------|-------|------|------|
| **知识** 📚 | cellwiki | `query_cell_type`, `list_cell_types`, `query_markers_for_genes` | 本地 CellWiki 知识库查询 |
| **聚类质量** 🔍 | clustering_quality | `silhouette_check`, `modularity_check`, `phiclust_check` | 聚类质量多工具评估 |
| **差异表达** 🧪 | clustering_degs | `find_markers`, `consensus_markers`, `resolution_scan` | 标记基因识别与分辨率扫描 |
| **格式转换** 🔄 | format_conversion | `load_h5ad`, `load_rds`, `load_mtx`, `load_sce` | 支持 h5ad/Seurat/mtx/SCE |
| **质量控制** 🛡️ | quality_control | `compute_qc_metrics`, `filter_cells`, `detect_doublets` | QC 指标、过滤、双细胞检测 |
| **批次校正** ⚖️ | batch_correction | `correct_harmony`, `correct_bbknn`, `correct_combat` | 三种主流批次校正方法 |
| **空间组学** 🗺️ | spatial | `load_visium`, `load_merfish`, `load_stereoseq`, `load_spatial` | 10x Visium / MERFISH / Stereo-seq |
| **轨迹推断** 🧭 | trajectory | `paga_trajectory` | PAGA 发育轨迹推断 |
| **可视化** 🎨 | plotting | `plot_embedding`, `violin_plot`, `dot_plot`, `heatmap` | 出版质量可视化（Okabe-Ito 色盲友好配色） |

### 4.2 Skills 设计特点

- **SKILL.md manifest** 📋: 每个 skill 目录含 YAML frontmatter，DeepAgents 自动发现
- **Lazy loading** ⚡: `skills/__init__.py` 使用 try/except 块按需加载，可选依赖（rpy2, harmonypy）缺失时优雅降级
- **纯函数模块** 🧩: MCTS 算法无 LLM/IO 依赖，便于测试和复用
- **SQLite 缓存** 💾: 证据缓存支持 TTL 过期，避免重复 LLM 查询

---

## 5. 技术亮点

### 5.1 MCTS 树搜索注释

传统方法：LLM 直接给出一个答案，可能置信度低。

Cell# 方法：
1. **Initialize**: 以 cluster 标记基因构建根节点假设 🌱
2. **Select**: UCB1 算法平衡探索与利用 🎯
3. **Expand**: 生成 2-3 个候选细胞类型假设 🌿
4. **Validate**: CellWiki 检索标记基因证据（Jaccard 相似度）🔍
5. **Evaluate**: 多维度评分（标记匹配 0.4 + 本体一致性 0.3 + 证据多样性 0.2 + LLM 置信度 0.1）📊
6. **Backpropagate**: 更新路径上所有节点统计 🔄
7. **Extract**: 返回最高平均奖励的节点作为最终标注 🏆

### 5.2 聚类质量多工具投票

单一指标容易误导（silhouette 倾向少聚类，modularity 倾向多聚类）。Cell# 采用多数投票制：
- **silhouette ASW > 0.5**: 接受；0.25-0.5: 中等；< 0.25: 调节 resolution
- **modularity Q > 0.3**: 接受；< 0.2: 调节 resolution
- **phiclust > 0 for >50% clusters**: 需要更高分辨率

当工具意见冲突时，采用多数投票制自动决策。🗳️

### 5.3 HITL 交互审核

低置信度标注（< 0.5）自动触发人工审核流程：
- 仅在全交互环境（TTY）下激活，notebook/CI 中自动跳过
- 用户可选择：✅ Accept / ❌ Reject / ✏️ Edit（自定义细胞类型名）
- 生成 `hitl_review.md` 审核摘要

### 5.4 配置系统四层优先级

1. **Defaults**: 代码中定义的默认值（Pydantic 模型）
2. **YAML 文件**: `config.yaml` 用户自定义配置
3. **环境变量**: `CELLHASHTAG_*` 前缀的环境变量
4. **运行时覆盖**: `CellHashtagAgent(profile="fast", llm__model="qwen-plus")`

支持 profile 切换：
- `fast` ⚡: max_iterations=1, max_anno_iter=2（快速验证）
- `default` 📦: 标准配置
- `deep` 🔬: 更多迭代和更深入的搜索

---

## 6. 输入输出

### 输入格式

| 格式 | 支持 | 说明 |
|------|------|------|
| `.h5ad` (AnnData) | ✅ 原生 | 推荐，最快 |
| `.rds` (Seurat) | ✅ rpy2 + zellkonverter | R Seurat 对象 |
| `.rds` (SCE) | ✅ zellkonverter | SingleCellExperiment 对象 |
| 10x `.mtx` | ✅ 原生 | matrix 目录 |
| 10x Visium | ✅ scanpy | 空间转录组 |
| MERFISH / Stereo-seq | ✅ CSV/parquet | 空间转录组 |

### 输出文件

| 文件 | 格式 | 说明 |
|------|------|------|
| `annotated.h5ad` | AnnData | 含 `Cell#` 和 `Cell#_confidence` 列 |
| `annotation_table.csv` | CSV | 每 cluster 一行：cluster, cell_type, confidence, reasoning |
| `annotation_report.md` | Markdown | 完整统计摘要 + 逐 cluster 表格 |
| `umap_celltype.png` | PNG | UMAP 按细胞类型着色 |
| `umap_confidence.png` | PNG | UMAP 按置信度着色（RdYlGn colormap） |
| `hitl_review.md` | Markdown | HITL 审核摘要（如有低置信度标注） |

---

## 7. 开发指标

| 指标 | 值 |
|------|-----|
| 总 Skills 数 | 8 (+ lats_search 测试) |
| 导出函数数 | 25 |
| 技术债状态 | ✅ 全部清零 (v3.1) |
| 测试覆盖 | MCTS 有测试，orchestrator/clustering/cellwiki 待补充 |
| 文档状态 | 📚 全面同步 |
| 支持数据格式 | 7 种 (h5ad, rds-Seurat, rds-SCE, 10x-mtx, Visium, MERFISH, Stereo-seq) |

---

## 8. 未来方向

- **Orchestrator Skills 集成** 🔗: 在 annotation 阶段自动调用新 skills（质量评估、DEG、批次校正）
- **并行化** ⚡: annotation 节点对 cluster 的串行调用 → 并行处理
- **测试覆盖** 🧪: 提升至 >70%
- **空间转录组分析** 🗺️: 在 annotation 中利用空间位置信息
- **自动 Skill Discovery** 🔎: 从 skills/ 目录自动发现，无需硬编码路径

---

## 附录：CellWiki 知识图谱项目

### 定位

CellWiki 是一个 **自维护的单细胞多组学知识库**，采用 "LLM Wiki" 模式——不是传统 RAG 在查询时从原始文档检索，而是用 LLM 增量构建和维护一个持久化的 Markdown wiki。知识在写入时编译一次，后续查询直接读取，无需每次重新推导。

> 🎯 **核心类比**：维基百科 vs 图书馆——CellWiki 是维基百科（知识已组织好），RAG 是图书馆（每次现找现翻）。

### 当前状态

| 指标 | 值 |
|------|-----|
| 版本 | v0.2 |
| 细胞类型 wiki 页面 | **93 个** |
| 总 wiki 内容 | 4,879 行 Markdown |
| 覆盖领域 | **肿瘤免疫学**（T cell exhaustion, tumor microenvironment） |
| 实体类型 | 6 种（cell_type, marker_gene, tissue, disease, method, trajectory） |

### 架构

```
Raw Sources (PDFs) ──► LLM Extraction ──► Wiki Pages (Markdown)
                           │                      │
                     Contradiction           Schema Rules
                     Detection               (schema.md)
                           │                      │
                     Knowledge Graph ◄─────── Multi-Omics Merger
```

### 四大 LangGraph 工作流

| Subgraph | 功能 | 关键步骤 |
|----------|------|----------|
| **Ingest** 📥 | 论文知识提取 | LLM 分析 → 人工审核 → wiki 生成 |
| **Query** 🔍 | 智能问答 | 实体匹配 → 图扩展 → 预算感知上下文 → LLM 综合 |
| **Lint** 🧹 | 自修复循环 | 扫描 → 分类 → 自动修复 → 复检 |
| **Research** 🔬 | 知识发现 | 生成查询 → 网络搜索 → 综合 → 提议更新 |

### 实体模型（6 种）

| 实体 | 说明 | 状态 |
|------|------|------|
| **cell_type** | Identity/state/context 三元组模型 | ✅ 93 页 |
| **marker_gene** | 特异性/敏感性/稳定性评分 | 规划中 |
| **tissue** | 细胞类型丰度映射 | 规划中 |
| **disease** | 细胞类型/基因关联 | 规划中 |
| **method** | 实验技术描述 | 规划中 |
| **trajectory** | 状态转换路径 | 规划中 |

### 证据层级（5 级）

| Tier | 定义 | 示例 |
|------|------|------|
| **1** | 直接实验验证（KO、功能实验、空间共定位） | FOXP3 KO 消除 Treg 功能 |
| **2** | 多组学一致性（RNA + 蛋白 + 表观一致） | scRNA-seq + CITE-seq + ATAC-seq 一致 |
| **3** | 单组学 + 多篇独立论文（≥3） | 3+ 篇报道相同标记 |
| **4** | 单篇论文报道 | 首次观察 |
| **5** | LLM 推理/假设/未验证 | 从共表达推断 |

### 矛盾处理机制

- **检测**：Ingest 时比较新提取与已有 wiki 数据
- **分类**：技术差异（不同平台/批次）vs 生物学差异（不同组织/疾病）
- **记录**：写入 `contradictions.md`，含来源、状态、涉及页面
- **消解**：技术→标注为平台差异；生物→细化实体定义；未知→触发 Deep Research

---

## CellWiki 与 Cell# 的关系

```
        CellWiki (知识库) ──供给──► Cell# (标注工具)
             📚                       🤖
        93 个细胞类型            利用 CellWiki 知识
        wiki 页面               做 MCTS 树搜索注释
```

- **CellWiki** 是 **知识供给端**：从论文中增量提取、整合、消解矛盾，维护结构化知识库
- **Cell#** 是 **知识消费端**：在 annotation 阶段通过 `cellwiki` skill 查询 CellWiki，辅助 LLM 做出专家级细胞类型判断
- **闭环**：Cell# 的新发现可以反哺 CellWiki（如新的标记基因、细胞亚型），CellWiki 的更新自动提升 Cell# 的注释质量
- **技术纽带**：Cell# 的 `skills/cellwiki/` 直接读取 CellWiki 的 `wiki/cell_types/*.md` 文件，通过 Jaccard 相似度匹配候选细胞类型

---

*Report generated: 2026-06-11 | Cell# v3.1 + CellWiki v0.2*

---

## 附录：Research Harness 系统

### 定位

`.research_harness_system/` 是一套 **渐进式披露（Progressive Disclosure）的 Agent 文档系统**，专为 LLM 辅助开发设计。核心思想：Agent 不应一次性读取所有文档，而是通过路由表按需加载，类似人类开发者的 "知道去哪找" 能力。

> 🎯 **核心类比**：Linux 的 `man` 命令——你不是读遍所有手册页，而是通过索引找到你需要的那一页。

### 架构

```
AGENTS.md (路由表)
    │
    ├── Level 1: 项目骨架（修改代码前必读）
    │       ARCHITECTURE.md, LEARNINGS.md, ACTIVE_VIBE.md, QUICKSTART.md
    │
    ├── Level 2: 领域知识（按任务类型加载）
    │       methods/*.md, workflows/*.md, SKILL_PLAN.md
    │
    ├── Level 3: 用户体验与输出（面向用户的变更）
    │       user-personas.md, output-format.md, hitl-spec.md, UX_SENSE.md
    │
    └── Level 4: 外部参考与 API 文档（查询框架/库 API 时加载）
            API 文档索引, references/*.md
```

### 核心文件（26 个文档）

| 文件 | 层级 | 说明 |
|------|------|------|
| **AGENTS.md** 🧭 | L0 | Agent 路由地图，唯一每次必读的文件 |
| **ARCHITECTURE.md** 🏗️ | L1 | v3.0 分层架构、模块边界、Breaking Changes |
| **LEARNINGS.md** 📝 | L1 | 已验证策略、反模式、知识盲区（只追加不覆盖） |
| **ACTIVE_VIBE.md** 📊 | L1 | 当前进度、待办清单、失败记录 |
| **QUICKSTART.md** 🚀 | L1 | 一页速览：文件位置、常用命令 |
| **TECHNICAL_DEBT.md** 💳 | L1 | 技术债登记表（按优先级追踪） |

### 领域知识文档（L2）

| 文档 | 内容 |
|------|------|
| **cell-biology-basis.md** 🧬 | 细胞类型层次、标记基因原则、多组学差异 |
| **clustering-quality.md** 🔍 | ASW/Modularity 决策逻辑 |
| **theoretical-basis.md** 📐 | 评分标准、分辨率调整公式 |
| **llm-prompt-principles.md** 🤖 | LLM Prompt 设计原则 |
| **scientific-beliefs.md** 🔬 | 多组学差异、注释不确定性处理 |
| **SKILL_PLAN.md** 📋 | v3.1 skills 路线图、函数签名 |

### 工作流与 UX 文档（L3）

| 文档 | 内容 |
|------|------|
| **langgraph-pipeline.md** | LangGraph State 字段、节点拓扑、数据流 |
| **hitl-spec.md** | 人机交互规范：低置信度注释的人工介入点 |
| **output-format.md** | AnnData 列定义、CSV/Markdown/UMAP 规范 |
| **user-personas.md** | 目标用户特征、错误处理要求 |
| **typical-analysis-pipeline.md** | 输入输出期望、上下游生态 |
| **empirical-defaults.md** | 数据处理、LLM、输出默认参数 |
| **UX_SENSE.md** | 进度显示、参数命名、静默模式 |

### 外部参考（L4）

| 文档 | 内容 |
|------|------|
| **API 文档索引** 📚 | 164 页提取的 LangGraph/LangChain/DeepAgent API |
| **references/index.md** | Scanpy、CellMarker、CellWiki 等外部参考 |
| **lats-framework.md** | MCTS 集成、Prompt 模板、配置 |
| **scanpy-llms.md** | AnnData 处理模式、序列化注意事项 |

### 决策表

Agent 根据任务类型自动决定加载哪些文档：

| 任务 | 必须加载 |
|------|----------|
| 修改 orchestrator | ARCHITECTURE + LangGraph 工作流 |
| 修改 clustering | ARCHITECTURE + clustering-quality |
| 修改 annotation | DeepAgents API + cell-biology-basis |
| 修改 prompt | llm-prompt-principles |
| 修复 bug | ACTIVE_VIBE（失败记录） + LEARNINGS |
| 首次介入 | QUICKSTART + ARCHITECTURE + ACTIVE_VIBE |

### 设计原则

1. **渐进式加载** 📦: 不一次性读取所有文档，按需加载
2. **只追加不覆盖** ✏️: LEARNINGS.md 和 TECHNICAL_DEBT.md 采用 Append-only 模式
3. **经验反哺** 🔄: 新发现追加到 LEARNINGS.md，缺陷追加到 TECHNICAL_DEBT.md
4. **路由表驱动** 🧭: AGENTS.md 是唯一入口，指引后续加载路径
5. **层级分离** 🔢: L0 路由 → L1 骨架 → L2 领域 → L3 UX → L4 外部 API

### 与 Cell# 的关系

```
Harness System (开发辅助) ──指导──► Cell# 开发流程
         🧭                           🤖
   26 个渐进式文档              LangGraph 标注 Agent
   按需加载                     4-node pipeline
   Agent 路由表                 8 个 skills, 25 个函数
```

Harness 系统确保开发 Agent 在正确的时间加载正确的上下文，避免知识缺失或知识过载。它是 Cell# 的 "开发基础设施"，而非运行时组件。

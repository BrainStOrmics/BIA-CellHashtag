"""
CellHashtagState — LangGraph 状态定义

v2.0 重构：
- AnnData 对象不放在 state 中，用文件路径代替
- 支持多组学、双层注释、多维度置信度
- 包含 clustering 迭代历史和质量评估
"""

from typing import TypedDict, Annotated, Optional
import operator


class CellHashtagState(TypedDict):
    """
    CellHashtag 的完整 LangGraph 状态。

    设计原则：
    - 大数据对象（AnnData）通过文件路径引用，不直接存入 state
    - 所有字段有明确的类型注解
    - anno_results 使用 operator.add 实现 map-reduce 聚合
    """

    # ===== 输入 =====
    input_path: str                     # 输入文件路径（h5ad / rds / mtx / csv）
    input_format: str                   # 输入格式：h5ad / rds / mtx / csv
    metadata_path: Optional[str]        # 可选：额外 metadata CSV 路径
    omics_type: str                     # scRNA / snRNA / spRNA / ATAC / CITE / Metab / multi
    cluster_key: str                    # adata.obs 中的聚类列名，默认 leiden

    # ===== LLM 配置 =====
    anno_llm: Optional[object]          # LangChain ChatOpenAI 实例
    max_iterations: int                 # 聚类迭代最大次数，默认 3
    max_anno_iter: int                  # 注释自批评最大迭代次数，默认 5

    # ===== Data Perception =====
    perceived_info: dict                # {species, platform, n_cells, n_genes, ...}
    adata_path: str                     # 标准化后的 AnnData pickle 路径

    # ===== Metadata =====
    metadata_dict: dict                 # {cluster_id: {metadata_key: value}}
    metadata_summary: str               # Metadata 文本摘要
    metadata_confirmed: bool            # HITL 是否已确认

    # ===== Clustering =====
    clustering_params: dict             # {resolution, n_neighbors, n_pcs, ...}
    clustering_result: dict             # {n_clusters, cluster_labels}
    clustering_history: list            # [{iteration, params, quality, decision}]
    clustering_quality: dict            # 质量评估结果 {phiclust, silhouette, modularity}
    iteration_count: int                # 当前聚类迭代次数

    # ===== Feature Extraction =====
    features: dict                      # {omics_type: {cluster_id: marker_list}}

    # ===== Knowledge Match =====
    matched_types: dict                 # {cluster_id: [(cell_type, score), ...]}
    contradictions: list                # 标记矛盾列表
    evidence_tiers: dict                # {cluster_id: evidence_tier (1-5)}

    # ===== Layer 1: Cell Type =====
    cell_type_annotations: dict         # {cluster_id: cell_type}
    layer1_confidence: dict             # {cluster_id: confidence_score (0-1)}

    # ===== Layer 2: Sub-annotation =====
    subtype_annotations: dict           # {cluster_id: [subtype_labels]}
    state_annotations: dict             # {cluster_id: [state_labels]}
    function_annotations: dict          # {cluster_id: [function_labels]}
    spatial_annotations: dict           # {cluster_id: [spatial_labels]}
    layer2_confidence: dict             # {cluster_id: {subtype: score, ...}}

    # ===== Confidence =====
    confidence: dict                    # 多维度置信度 {overall, marker_match, ...}

    # ===== Annotation Sub-Graph (per cluster) =====
    # 这些字段用于 subgraph 的 per-cluster 处理
    current_cluster: str                # 当前正在处理的 cluster ID
    anno_strategy: str                  # primary / secondary
    cell_markers_table: str             # Markdown 表格：cell types 和 markers
    exp_summary_dict: dict              # {cluster_id: expression_summary_markdown}
    anno_result: Annotated[list, operator.add]  # 单个 cluster 的注释结果
    anno_celltype: str                  # 当前 cluster 的注释结果
    anno_critique: str                  # 自批评反馈
    critic_decision: str                # Approved / Disapproved
    n_anno_iter: int                    # 当前注释迭代次数

    # ===== Output =====
    hashtags: dict                      # {cluster_id: ["#CD4_T_cell", "#regulatory", ...]}
    report_path: str                    # Markdown 报告输出路径
    annotation_table_path: str          # CSV 表格输出路径
    output_adata_path: str              # 输出 AnnData pickle 路径

    # ===== Status =====
    status: str                         # perception → metadata → clustering → ... → done
    errors: list                        # 错误消息列表
    interrupt_message: str              # HITL 交互消息

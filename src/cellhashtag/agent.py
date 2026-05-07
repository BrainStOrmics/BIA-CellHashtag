"""
CellHashtag v2.0 — 主 Agent 类

重构自 _Agent.py，主要变化：
- 内置 Clustering + 质量评估迭代循环
- 多格式输入支持
- 双层注释体系
- CellWiki 知识集成（stub）
- HITL 交互支持
- 丰富输出（AnnData + CSV + Markdown 报告）
"""

import os
import pickle
from typing import Optional, Any
from pathlib import Path

# LangGraph
from langgraph.checkpoint.memory import MemorySaver
from langgraph.pregel import RetryPolicy

# 内部模块
from .state import CellHashtagState
from .config.config import LLMConfig, AgentConfig, load_config, setup_llm
from .utils.cellwiki_client import CellWikiClient


class CellHashtagAgent:
    """
    CellHashtag v2.0 Agent。

    基于 LangGraph 的单细胞注释 Agent，支持多组学、
    内置聚类质量评估和双层注释。
    """

    def __init__(
        self,
        config_yaml: str = None,
        web_scraper_api_key: Optional[str] = None,
        cellwiki_dir: Optional[str] = None,
    ):
        """
        初始化 CellHashtag Agent。

        Args:
            config_yaml: 配置文件路径。默认使用 config/config.yaml。
            web_scraper_api_key: Web 搜索 API key（可选）。
            cellwiki_dir: CellWiki wiki/ 目录路径（可选）。
        """
        print("Initializing Cell# v2.0 Agent")

        # 加载配置
        if config_yaml:
            load_config(config_yaml)
        else:
            default_yaml = str(Path(__file__).parent / "config" / "config.yaml")
            if os.path.exists(default_yaml):
                load_config(default_yaml)

        self.web_scraper_api_key = web_scraper_api_key
        self.cellwiki = CellWikiClient(wiki_dir=cellwiki_dir)
        self.graph = None
        self.memory = MemorySaver()

        # 初始化图
        try:
            self.graph = self._build_graph()
            print("Cell# v2.0 agent ready.")
        except Exception as e:
            print(f"Failed to initialize graph: {e}")
            raise

    def _build_graph(self):
        """构建 LangGraph 图。"""
        from .graphs.orchestrator import build_orchestrator_graph

        builder = build_orchestrator_graph()
        return builder.compile(checkpointer=self.memory)

    def run(
        self,
        input_path: str,
        llm: Any = None,
        input_format: str = "h5ad",
        cluster_key: str = "leiden",
        omics_type: str = "auto",
        metadata_path: str = None,
        cell_marker_df_dir: str = "",
        max_iterations: int = 3,
        max_anno_iter: int = 5,
        thread_id: str = "1",
        output_dir: str = "output",
        plot_annotation: bool = True,
    ) -> dict:
        """
        运行 CellHashtag 注释流程。

        Args:
            input_path: 输入文件路径（.h5ad / .rds / mtx 目录）。
            llm: LangChain ChatOpenAI 实例。如果为 None，从配置初始化。
            input_format: 输入格式：h5ad / rds / mtx。
            cluster_key: adata.obs 中的聚类列名。
            omics_type: 组学类型，"auto" 表示自动推断。
            metadata_path: 额外 metadata CSV 路径（可选）。
            cell_marker_df_dir: 细胞标记基因 CSV 路径。
            max_iterations: 聚类最大迭代次数。
            max_anno_iter: 注释自批评最大迭代次数。
            thread_id: LangGraph thread ID。
            output_dir: 输出目录。
            plot_annotation: 是否绘制 UMAP 注释图。

        Returns:
            完整的结果字典。
        """
        # 初始化 LLM
        if llm is None:
            llm = setup_llm()
            if llm is None:
                raise ValueError(
                    "LLM not configured. Provide llm parameter or set up config.yaml."
                )

        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)

        # 加载细胞标记基因
        cell_markers_table = self._load_cell_markers(cell_marker_df_dir)

        # 构建初始状态
        initial_state = {
            "input_path": input_path,
            "input_format": input_format,
            "metadata_path": metadata_path,
            "omics_type": omics_type,
            "cluster_key": cluster_key,
            "anno_llm": llm,
            "max_iterations": max_iterations,
            "max_anno_iter": max_anno_iter,
            "cell_markers_table": cell_markers_table,
            "adata_path": os.path.join(output_dir, "adata_temp.pkl"),
            "output_adata_path": os.path.join(output_dir, "annotated_adata.pkl"),
            "report_path": os.path.join(output_dir, AgentConfig.REPORT_FILENAME),
            "annotation_table_path": os.path.join(output_dir, AgentConfig.TABLE_FILENAME),
            # 初始化空字段
            "perceived_info": {},
            "metadata_dict": {},
            "metadata_summary": "",
            "metadata_confirmed": False,
            "clustering_params": {"resolution": 0.5},
            "clustering_result": {},
            "clustering_history": [],
            "clustering_quality": {},
            "iteration_count": 0,
            "features": {},
            "matched_types": {},
            "contradictions": [],
            "evidence_tiers": {},
            "cell_type_annotations": {},
            "layer1_confidence": {},
            "subtype_annotations": {},
            "state_annotations": {},
            "function_annotations": {},
            "spatial_annotations": {},
            "layer2_confidence": {},
            "confidence": {},
            "hashtags": {},
            "status": "initializing",
            "errors": [],
            "interrupt_message": "",
            "current_cluster": "",
            "anno_strategy": "primary",
            "exp_summary_dict": {},
            "anno_result": [],
            "anno_celltype": "",
            "anno_critique": "",
            "critic_decision": "",
            "n_anno_iter": 0,
        }

        config = {"configurable": {"thread_id": thread_id}}

        # 执行图
        print("Starting Cell# annotation pipeline...")
        result = self.graph.invoke(input=initial_state, config=config)

        # 后处理输出
        result = self._post_process(result, output_dir, plot_annotation)

        return result

    def _load_cell_markers(self, marker_path: str) -> str:
        """加载细胞标记基因为 Markdown 表格。"""
        import pandas as pd

        if not marker_path or not os.path.exists(marker_path):
            return "No cell markers provided."

        try:
            df = pd.read_csv(marker_path)
            # 假设格式: cell_type, cell_markers
            lines = ["| Cell Type | Markers |", "|-----------|---------|"]
            for _, row in df.iterrows():
                cell_type = row.get("cell_type", "Unknown")
                markers = row.get("cell_markers", "")
                lines.append(f"| {cell_type} | {markers} |")
            return "\n".join(lines)
        except Exception as e:
            print(f"Warning: Failed to load cell markers: {e}")
            return "Failed to load cell markers."

    def _post_process(self, result: dict, output_dir: str, plot: bool) -> dict:
        """后处理：生成报告、表格、可视化。"""
        try:
            import pickle
            with open(result.get("adata_path", ""), "rb") as f:
                adata = pickle.load(f)

            # 更新 AnnData.obs
            cell_types = result.get("cell_type_annotations", {})
            cluster_key = result.get("cluster_key", "leiden")

            if cell_types:
                ct_map = {}
                for cluster, ct in cell_types.items():
                    ct_map[cluster] = ct

                adata.obs["Cell#"] = adata.obs[cluster_key].map(ct_map)
                adata.obs["Cell#_raw"] = adata.obs[cluster_key].map(ct_map)

                # 保存
                with open(result.get("output_adata_path", ""), "wb") as f:
                    pickle.dump(adata, f)

            # 绘图
            if plot:
                self._plot_annotation(adata, cluster_key)

        except Exception as e:
            print(f"Warning: Post-processing failed: {e}")

        return result

    def _plot_annotation(self, adata, cluster_key: str):
        """绘制 UMAP 注释图。"""
        try:
            import scanpy as sc
            if "X_umap" not in adata.obsm:
                print("Computing UMAP...")
                sc.pp.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)

            sc.pl.umap(
                adata,
                color=[cluster_key, "Cell#", "Cell#_raw"],
                ncols=1,
                save="_cellhashtag.pdf",
            )
        except Exception as e:
            print(f"Warning: Failed to plot: {e}")

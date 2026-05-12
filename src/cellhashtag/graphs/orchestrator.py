"""
CellHashtag 主图编排器。

构建完整的 LangGraph 工作流：
  Data Perception → Metadata → Clustering → Feature Extraction →
  Knowledge Match → Annotation L1 → Annotation L2 → Confidence → Output
"""

import os
import json
import re
import pickle
import importlib
from typing import Any

from langgraph.graph import StateGraph, START, END
from langgraph.types import Send
from langchain_core.prompts import PromptTemplate
from langchain_openai import ChatOpenAI

from ..state import CellHashtagState
from ..utils.utilities import get_exp_summary


def build_orchestrator_graph() -> StateGraph:
    """
    构建主 LangGraph 编排器。

    Returns:
        编译后的 LangGraph。
    """
    builder = StateGraph(CellHashtagState)

    # 注册节点
    builder.add_node("data_perception", _import_node("data_perception", "node_data_perception"))
    builder.add_node("metadata_extraction", _import_node("metadata_node", "node_metadata_extraction"))
    builder.add_node("clustering", _import_node("clustering_node", "node_clustering"))

    # 边：线性流程
    builder.add_edge(START, "data_perception")
    builder.add_edge("data_perception", "metadata_extraction")
    builder.add_edge("metadata_extraction", "clustering")

    # 条件边：聚类可能需要重试
    builder.add_conditional_edges(
        "clustering",
        _clustering_router,
        {
            "clustering_retry": "clustering",
            "clustering_done": "annotation_subgraph",
        },
    )

    # Annotation subgraph（map-reduce）
    builder.add_node("annotation_subgraph", _build_annotation_subgraph_node())

    builder.add_edge("annotation_subgraph", END)

    return builder


def _import_node(module_name: str, func_name: str):
    """动态导入节点函数。"""
    module = importlib.import_module(f"cellhashtag.nodes.{module_name}")
    return getattr(module, func_name)


def _clustering_router(state: dict) -> str:
    """聚类路由决策。"""
    status = state.get("status", "")
    if status == "clustering_retry":
        return "clustering_retry"
    return "clustering_done"


def _build_annotation_subgraph_node():
    """
    Build annotation subgraph node using LATS tree search.
    
    For each cluster, delegates to low_hierarchy graph which uses
    LATS as the default search strategy.
    """
    
    def annotation_subgraph_node(state: dict) -> dict:
        """Map-Reduce annotation with LATS search."""
        # Load adata
        with open(state["adata_path"], "rb") as f:
            adata = pickle.load(f)
        
        cluster_key = state.get("cluster_key", "leiden")
        clusters = adata.obs[cluster_key].unique().tolist()
        
        # Load config for LATS
        config_path = os.path.join(os.path.dirname(__file__), "..", "config", "config.yaml")
        config = load_config(config_path) if os.path.exists(config_path) else {}
        
        # Build low hierarchy graph (LATS-based)
        low_graph = create_low_hierarchy_graph(config)
        compiled_low = low_graph.compile()
        
        all_results = []
        
        for cluster in clusters:
            # Prepare LATS state for this cluster
            lats_input = {
                "cluster_id": str(cluster),
                "cluster_markers": _get_cluster_markers(adata, cluster_key, cluster),
                "cluster_expression_summary": _get_cluster_expression(adata, [str(cluster)]),
                "tissue_source": state.get("perceived_info", {}).get("tissue", "unknown"),
                "experimental_condition": state.get("metadata_summary", ""),
                
                # Constraints from harness
                "allowed_cell_types": None,  # Could filter based on tissue
                "excluded_markers": [],
                "priority_pathways": [],
                
                # Config paths
                "lats_config_path": config_path,
                "lats_prompts_dir": os.path.join(os.path.dirname(__file__), "..", "prompts", "low_hierarchy"),
                
                # Iteration control
                "annotation_iteration": 0,
                "max_annotation_iter": state.get("max_anno_iter", 5),
                "min_confidence": config.get("annotation", {}).get("confidence_threshold", 0.7),
            }
            
            # Execute LATS search
            result = compiled_low.invoke(lats_input)
            
            # Extract annotation result
            annotation = result.get("annotation_result")
            if annotation and not result.get("fallback_to_manual"):
                all_results.append({
                    "cluster": str(cluster),
                    "cell_type": annotation.get("cell_type", "Unknown"),
                    "reason": "; ".join(annotation.get("reasoning", [])),
                    "confidence": annotation.get("confidence", 0.5),
                    "evidence": annotation.get("evidence", []),
                })
            else:
                # Fallback to simple annotation if LATS fails
                fallback = _simple_annotation_fallback(cluster, adata, cluster_key, state)
                all_results.append(fallback)
        
        # Aggregate results
        cell_type_annotations = {r["cluster"]: r["cell_type"] for r in all_results}
        layer1_confidence = {r["cluster"]: r["confidence"] for r in all_results}
        
        return {
            "cell_type_annotations": cell_type_annotations,
            "layer1_confidence": layer1_confidence,
            "hashtags": {cid: [f"#{ct}"] for cid, ct in cell_type_annotations.items()},
            "annotation_details": {r["cluster"]: r for r in all_results},
            "status": "annotation_done",
        }
    
    return annotation_subgraph_node


def _get_cluster_markers(adata, cluster_key: str, cluster) -> list:
    """Extract top marker genes for a cluster."""
    try:
        import scanpy as sc
        if "rank_genes_groups" not in adata.uns:
            sc.tl.rank_genes_groups(adata, groupby=cluster_key, reference="rest")
        markers = list(adata.uns["rank_genes_groups"]["names"][cluster][:20])
        return [m for m in markers if m and not m.startswith("MT-")]
    except Exception:
        return []


def _get_cluster_expression(adata, clusters) -> Dict[str, float]:
    """Get average expression summary for clusters."""
    try:
        return get_exp_summary(adata, clusters)
    except Exception:
        return {}


def _simple_annotation_fallback(cluster, adata, cluster_key, state) -> dict:
    """Fallback annotation when LATS search doesn't produce confident result."""
    # Simple marker-based annotation as fallback
    try:
        import scanpy as sc
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, reference="rest")
        top_genes = list(adata.uns["rank_genes_groups"]["names"][cluster][:10])
        
        # Basic heuristic: match against known markers
        cell_markers = state.get("cell_markers_table", "")
        for line in cell_markers.split("
")[2:]:  # Skip header
            if "|" in line:
                parts = line.split("|")
                if len(parts) >= 3:
                    ct, markers = parts[1].strip(), parts[2].strip()
                    if any(m.strip() in top_genes for m in markers.split(",")):
                        return {
                            "cluster": str(cluster),
                            "cell_type": ct,
                            "reason": f"Marker match: {top_genes[:5]}",
                            "confidence": 0.5,
                            "evidence": [],
                        }
    except Exception:
        pass
    
    return {
        "cluster": str(cluster),
        "cell_type": "Unknown",
        "reason": "Fallback: no confident annotation found",
        "confidence": 0.3,
        "evidence": [],
    }




def _annotate_cluster(cluster, adata, cluster_key, llm, prompt, critic_prompt, state):
    """单个 cluster 的 annotation → critic 循环。"""
    # 提取 cluster 的 marker 基因
    try:
        import scanpy as sc
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, reference="rest")
        top_genes = list(adata.uns["rank_genes_groups"]["names"][cluster][:15])
    except Exception:
        top_genes = []

    # 生成表达摘要
    exp_summary = _get_cluster_expression(adata, top_genes)
    metadata = state.get("metadata_dict", {}).get(str(cluster), {})
    metadata_table = "\n".join(f"- {k}: {v}" for k, v in metadata.items())
    cell_markers = state.get("cell_markers_table", "No markers provided.")

    max_iter = state.get("max_anno_iter", 5)
    cell_type = "Unknown"
    reason = ""

    for i in range(max_iter):
        critique = ""
        if i > 0:
            critique = "\nPrevious critique: " + state.get("anno_critique", "")

        try:
            prompt_obj = PromptTemplate(
                input_variables=["exp_summary_table", "metadata_table", "cell_markers_table", "critique"],
                template=prompt,
            )
            chain = prompt_obj | llm
            response = chain.invoke({
                "exp_summary_table": exp_summary,
                "metadata_table": metadata_table,
                "cell_markers_table": cell_markers,
                "critique": critique,
            })

            cell_type, reason = _parse_annotation_response(response.content)

            # 自批评
            critic_obj = PromptTemplate(
                input_variables=["cell_type", "exp_summary_table", "metadata_table",
                                 "draft_response_wreason", "cell_markers_table"],
                template=critic_prompt,
            )
            critic_chain = critic_obj | llm
            critic_response = critic_chain.invoke({
                "cell_type": cell_type,
                "exp_summary_table": exp_summary,
                "metadata_table": metadata_table,
                "draft_response_wreason": reason,
                "cell_markers_table": cell_markers,
            })

            critic_decision = _parse_critic_response(critic_response.content)
            state["anno_critique"] = critic_response.content

            if critic_decision.lower() == "approved":
                break

        except Exception as e:
            print(f"Annotation error for cluster {cluster}: {e}")
            break

    return {
        "cluster": str(cluster),
        "cell_type": cell_type,
        "reason": reason,
        "confidence": 0.7 if cell_type != "Unknown" else 0.3,
    }


def _get_cluster_expression(adata, genes):
    """生成 cluster 的表达摘要。"""
    try:
        return get_exp_summary(adata, genes)
    except Exception:
        sep = ", "
        return f"Genes: {sep.join(genes[:10])}"


def _parse_annotation_response(content: str) -> tuple:
    """解析 LLM 的注释响应。"""
    json_match = re.search(r"```json\s*\[(.*?)\]\s*```", content, re.DOTALL)
    if json_match:
        try:
            result = json.loads(f"[{json_match.group(1)}]")
            return result[0], content[:200]
        except json.JSONDecodeError:
            pass
    return "Unknown", content[:200]


def _parse_critic_response(content: str) -> str:
    """解析 critic 响应。"""
    json_match = re.search(r"```json\s*\{(.*?)\}\s*```", content, re.DOTALL)
    if json_match:
        try:
            result = json.loads(f"{{{json_match.group(1)}}}")
            return result.get("final decision", "Disapproved")
        except json.JSONDecodeError:
            pass
    return "Disapproved"


def _load_prompt(prompt_name: str) -> str:
    """加载 prompt 模板。"""
    prompt_path = os.path.join(os.path.dirname(__file__), "..", "prompts", f"{prompt_name}.md")
    try:
        with open(prompt_path, "r", encoding="utf-8") as f:
            content = f.read()
        return content.replace("<<", "{").replace(">>", "}")
    except FileNotFoundError:
        return ""

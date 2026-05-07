"""
Metadata 提取节点。

从 AnnData.obs 和可选的额外 metadata CSV 中提取与聚类相关的 Metadata。
为 LLM 注释提供上下文信息。
"""

import pandas as pd


def node_metadata_extraction(state: dict) -> dict:
    """
    提取 Metadata。

    Args:
        state: CellHashtagState。

    Returns:
        更新后的 state 片段。
    """
    import pickle

    # 加载 AnnData
    with open(state["adata_path"], "rb") as f:
        adata = pickle.load(f)

    cluster_key = state.get("cluster_key", "leiden")
    clusters = adata.obs[cluster_key].unique().tolist()

    # 从 adata.obs 提取 metadata
    metadata_dict = {}
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        cluster_obs = adata.obs[mask]

        meta = {}
        for col in cluster_obs.columns:
            if col == cluster_key:
                continue
            # 数值列：计算均值和标准差
            if pd.api.types.is_numeric_dtype(cluster_obs[col]):
                mean_val = cluster_obs[col].mean()
                meta[f"{col}_mean"] = round(float(mean_val), 2)
            # 分类列：计算主要类别
            else:
                top_val = cluster_obs[col].value_counts().index[0]
                meta[f"{col}_top"] = str(top_val)

        metadata_dict[str(cluster)] = meta

    # 如果有额外 metadata CSV，合并
    metadata_path = state.get("metadata_path")
    if metadata_path:
        metadata_dict = _merge_external_metadata(metadata_dict, metadata_path, cluster_key)

    # 生成摘要
    metadata_summary = _generate_metadata_summary(metadata_dict)

    return {
        "metadata_dict": metadata_dict,
        "metadata_summary": metadata_summary,
        "metadata_confirmed": False,  # 等待 HITL 确认
        "status": "metadata_done",
    }


def _merge_external_metadata(
    metadata_dict: dict,
    metadata_path: str,
    cluster_key: str,
) -> dict:
    """合并外部 metadata CSV。"""
    try:
        ext_meta = pd.read_csv(metadata_path)
        # 假设 CSV 有 cluster_id 列
        if "cluster_id" in ext_meta.columns:
            for _, row in ext_meta.iterrows():
                cid = str(row["cluster_id"])
                if cid in metadata_dict:
                    for col in ext_meta.columns:
                        if col != "cluster_id":
                            metadata_dict[cid][col] = str(row[col])
    except Exception as e:
        print(f"Warning: Failed to load external metadata: {e}")

    return metadata_dict


def _generate_metadata_summary(metadata_dict: dict) -> str:
    """生成 Metadata 文本摘要。"""
    lines = ["## Metadata Summary\n"]
    for cluster_id, meta in metadata_dict.items():
        lines.append(f"### Cluster {cluster_id}")
        for key, val in meta.items():
            lines.append(f"- {key}: {val}")
        lines.append("")
    return "\n".join(lines)

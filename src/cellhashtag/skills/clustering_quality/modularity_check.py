"""
Modularity 模块度检查。

使用 igraph 计算 Leiden/Louvain 聚类的模块度。
模块度 ∈ [0, 1]:
  > 0.3  → 有意义的社区结构
  > 0.5  → 强社区结构
"""

from typing import Any


def modularity_check(
    adata: Any,
    cluster_key: str = "leiden",
) -> dict:
    """
    计算聚类模块度。

    Args:
        adata: AnnData 对象（需要 .obsp["connectivities"]）。
        cluster_key: adata.obs 中的聚类列名。

    Returns:
        {
            "tool": "modularity",
            "modularity": float,
            "recommendation": str,
            "details": str,
        }
    """
    try:
        import igraph as ig
    except ImportError:
        raise ImportError(
            "igraph is required. Install with: pip install python-igraph"
        )

    if "connectivities" not in adata.obsp:
        raise ValueError(
            "connectivities not found in adata.obsp. "
            "Run sc.pp.neighbors(adata) first."
        )

    # 构建 igraph 图
    conn = adata.obsp["connectivities"]
    rows, cols = conn.nonzero()
    weights = conn[rows, cols].A1 if hasattr(conn[rows, cols], "A1") else conn[rows, cols]

    g = ig.Graph(n=adata.n_obs, edges=list(zip(rows, cols)), edge_attrs={"weight": weights})

    # 计算模块度
    membership = adata.obs[cluster_key].astype("category").cat.codes.tolist()
    mod = g.modularity(membership, weights="weight")

    # 决策
    if mod < 0.2:
        recommendation = "adjust"
    elif mod > 0.6:
        recommendation = "ok"
    else:
        recommendation = "ok"

    return {
        "tool": "modularity",
        "modularity": round(float(mod), 4),
        "recommendation": recommendation,
        "details": f"Modularity={mod:.3f}",
    }

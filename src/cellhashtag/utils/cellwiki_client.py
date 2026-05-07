"""
CellWiki 查询接口。

CellHashtag 通过此接口查询 CellWiki 知识库，
获取标记基因、细胞类型、轨迹等信息用于注释。

当前为 stub 实现，后续与 CellWiki 项目对接。
"""

import os
from typing import Optional


class CellWikiClient:
    """CellWiki 知识查询客户端。"""

    def __init__(self, wiki_dir: str = None):
        """
        初始化 CellWiki 客户端。

        Args:
            wiki_dir: CellWiki wiki/ 目录路径。
        """
        if wiki_dir is None:
            # 默认路径：~/CellWiki/wiki/
            wiki_dir = os.path.expanduser("~/CellWiki/wiki")
        self.wiki_dir = wiki_dir
        self._cache = {}

    def query_marker_genes(self, gene_names: list = None) -> list:
        """
        查询标记基因信息。

        Args:
            gene_names: 基因名列表。如果为 None，返回所有。

        Returns:
            标记基因信息列表。
        """
        # TODO: 实际查询 CellWiki 的 marker_genes 页面
        return []

    def query_cell_types(self, cell_type: str = None) -> list:
        """
        查询细胞类型信息。

        Args:
            cell_type: 细胞类型名称。如果为 None，返回所有。

        Returns:
            细胞类型信息列表。
        """
        # TODO: 实际查询 CellWiki 的 cell_types 页面
        return []

    def match_cell_type(self, markers: list, omics_type: str = "scRNA") -> list:
        """
        根据标记基因匹配最可能的细胞类型。

        Args:
            markers: 观察到的标记基因列表。
            omics_type: 组学类型。

        Returns:
            [(cell_type, match_score), ...] 按分数降序排列。
        """
        # TODO: 实现匹配算法
        return []

    def get_subtype_candidates(self, cell_type: str, metadata: dict) -> list:
        """
        根据细胞类型和 Metadata 获取亚群候选。

        Args:
            cell_type: 基础细胞类型。
            metadata: 样本 Metadata。

        Returns:
            亚群候选列表。
        """
        # TODO: 查询 subtypes_by_context
        return []

    def check_contradictions(self, markers: list) -> list:
        """
        检查标记基因之间是否有已知矛盾。

        Args:
            markers: 标记基因列表。

        Returns:
            矛盾列表。
        """
        # TODO: 查询 contradictions.md
        return []

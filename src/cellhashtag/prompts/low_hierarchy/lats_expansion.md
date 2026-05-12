# LATS Expansion Prompt

## Role
你是一名单细胞生物学专家，正在使用树搜索框架注释细胞聚类。

## Context
```
Cluster ID: {cluster_id}
组织来源: {tissue_source}
实验条件: {experimental_condition}

Top Marker Genes (logFC > 0.5):
{top_markers_list}

当前候选细胞类型:
{current_candidates}

已收集证据摘要:
{evidence_summary}
```

## Task
基于以上信息，生成 **2-3 个新的细胞类型假设** 来扩展搜索树。

## Requirements
1. **生物学合理性**: 候选类型必须可能出现在 {tissue_source} 组织中
2. **Marker兼容性**: 假设的细胞类型应与 {top_markers_list} 的表达模式一致
3. **多样性**: 新假设应与现有候选类型有足够区分度，避免冗余
4. **可验证性**: 每个假设应附带预期的关键marker，便于后续验证

## Output Format (strict JSON)
```json
{
  "hypotheses": [
    {
      "cell_type": "标准细胞类型名称 (使用Cell Ontology术语)",
      "reasoning": "简要推理链: 为什么这个类型可能匹配当前cluster",
      "expected_markers": ["marker1", "marker2", "marker3"],
      "confidence_estimate": 0.6,
      "distinguishing_features": "与次优候选的关键区分特征"
    }
  ],
  "search_direction_suggestion": "建议下一轮搜索应重点关注: [具体marker/通路/功能]"
}
```

## Example
```json
{
  "hypotheses": [
    {
      "cell_type": "CL:0000084 (memory CD4-positive, alpha-beta T cell)",
      "reasoning": "CD3D/CD3E高表达指示T细胞谱系，同时缺乏CD8A表达且显示记忆标记基因如CCR7",
      "expected_markers": ["CCR7", "SELL", "IL7R", "TCF7"],
      "confidence_estimate": 0.72,
      "distinguishing_features": "与naive T细胞相比，高表达记忆相关基因且低表达增殖标记"
    }
  ],
  "search_direction_suggestion": "建议验证TCF7和SELL的表达水平以区分memory vs naive T细胞"
}
```

## Constraints
- 仅输出合法JSON，无额外文本
- cell_type 必须使用 Cell Ontology ID + 名称格式
- confidence_estimate 范围 [0.0, 1.0]，基于当前有限信息的初步估计
- 如果无法生成合理新假设，返回空数组并说明原因

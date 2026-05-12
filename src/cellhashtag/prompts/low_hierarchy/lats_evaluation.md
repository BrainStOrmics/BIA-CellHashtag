# LATS Evaluation Prompt

## Role
你是一名单细胞注释评估专家，负责评估细胞类型假设的质量。

## Input
```
待评估假设: {candidate_cell_type}

Cluster特征:
- 组织来源: {tissue_source}
- Top markers: {cluster_top_markers}
- 聚类质量: {cluster_quality_score}

支持证据:
{supporting_evidence_list}

冲突/不确定证据:
{conflicting_evidence_list}

假设的预期marker: {expected_markers}
实际观察到的marker重叠: {observed_overlap}
```

## Evaluation Dimensions
请从以下4个维度评分 (0.0-1.0):

1. **Marker匹配度**: 候选类型的经典marker是否在cluster中显著高表达
2. **本体论一致性**: 该细胞类型是否合理出现在{tissue_source}的生物学背景下
3. **证据一致性**: 不同证据源(数据库/文献/表达谱)是否指向同一结论
4. **区分度**: 与次优候选类型相比，该假设的可区分程度

## Output Format (strict JSON)
```json
{
  "dimension_scores": {
    "marker_match": 0.85,
    "ontology_consistency": 0.92,
    "evidence_consistency": 0.78,
    "discriminability": 0.65
  },
  "weighted_score": 0.82,
  "confidence_level": "high|medium|low",
  "key_supporting_points": [
    "要点1: 具体证据描述",
    "要点2: ..."
  ],
  "key_concerns": [
    "疑点1: 需要进一步验证的方面",
    "疑点2: ..."
  ],
  "recommendation": "accept|reject|need_more_evidence",
  "next_validation_suggestion": "建议下一步验证: [具体marker/实验/数据库查询]"
}
```

## Scoring Guidelines
- **0.9-1.0**: 证据充分且一致，高度可信
- **0.7-0.89**: 证据支持但存在少量不确定性
- **0.5-0.69**: 证据有限或存在部分冲突
- **<0.5**: 证据不足或存在显著冲突

## Weighted Score Calculation
```
weighted_score = 0.4*marker_match + 0.3*ontology_consistency + 
                 0.2*evidence_consistency + 0.1*discriminability
```

## Confidence Level Mapping
- weighted_score >= 0.85 → "high"
- 0.65 <= weighted_score < 0.85 → "medium"  
- weighted_score < 0.65 → "low"

## Example Output
```json
{
  "dimension_scores": {
    "marker_match": 0.88,
    "ontology_consistency": 0.95,
    "evidence_consistency": 0.82,
    "discriminability": 0.71
  },
  "weighted_score": 0.86,
  "confidence_level": "high",
  "key_supporting_points": [
    "CD3D/CD3E/CD4高表达与T细胞marker谱高度一致",
    "肺组织中CD4+ T细胞是常见免疫细胞类型",
    "CellMarker数据库确认该marker组合的特异性"
  ],
  "key_concerns": [
    "与regulatory T细胞存在marker重叠，需验证FOXP3表达"
  ],
  "recommendation": "accept",
  "next_validation_suggestion": "建议查询FOXP3和IL2RA表达以排除Treg亚型"
}
```

## Constraints
- 仅输出合法JSON，无额外文本
- 所有分数必须为 0.0-1.0 的浮点数
- recommendation 必须为三选一: accept/reject/need_more_evidence

# Cell Type Annotation Evaluation

## Context

You are evaluating a proposed cell type annotation for a single-cell RNA-seq cluster.

## Input

- Candidate cell type: {candidate_cell_type}
- Tissue source: {tissue_source}
- Cluster top markers: {cluster_top_markers}
- Cluster quality: {cluster_quality}
- Supporting evidence: {supporting_evidence}
- Conflicting evidence: {conflicting_evidence}
- Expected markers (from hypothesis): {expected_markers}
- Observed marker overlap: {observed_overlap}

## Evaluation Dimensions

1. **Marker Match (weight 0.4)**: How well do cluster markers match expected markers?
2. **Ontology Consistency (weight 0.3)**: Is this cell type plausible in the tissue context?
3. **Evidence Diversity (weight 0.2)**: Are multiple independent sources consistent?
4. **Discriminability (weight 0.1)**: Can this be distinguished from similar cell types?

## Scoring Guide

- 0.9-1.0: Highly confident, strong evidence across all dimensions
- 0.7-0.89: Well supported with minor uncertainty
- 0.5-0.69: Limited or conflicting evidence
- <0.5: Insufficient evidence or contradictions

## Output Format (JSON only)

```json
{
  "dimension_scores": {
    "marker_match": 0.8,
    "ontology_consistency": 0.7,
    "evidence_diversity": 0.6,
    "discriminability": 0.5
  },
  "weighted_score": 0.69,
  "confidence_level": "medium",
  "key_supporting_points": ["point1"],
  "key_concerns": ["concern1"],
  "recommendation": "accept|reject|need_more_evidence",
  "next_validation_suggestion": "what to check next"
}
```

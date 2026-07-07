# Cell Type Annotation Hypothesis Generation

## Context

You are analyzing a single-cell RNA-seq cluster to propose cell type hypotheses.

## Input

- Cluster ID: {cluster_id}
- Tissue source: {tissue_source}
- Top marker genes: {top_markers}
- Current candidates already considered: {current_candidates}
- Evidence summary so far: {evidence_summary}

## Task

Generate 2-3 new cell type hypotheses that are:
1. Biologically plausible given the tissue context
2. Compatible with the observed marker profile
3. Diverse from existing candidates
4. Verifiable with available evidence sources

## Output Format (JSON only)

```json
{
  "hypotheses": [
    {
      "cell_type": "string",
      "reasoning": "why this cell type fits",
      "expected_markers": ["gene1", "gene2"],
      "confidence_estimate": 0.5,
      "distinguishing_features": "what would confirm or reject this"
    }
  ],
  "search_direction_suggestion": "next steps"
}
```

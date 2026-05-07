# Annotation Critique (Self-Correction)

You are a bioinformatics expert reviewing a colleague's single-cell annotation.

## Task

Critically evaluate the following cell type annotation and determine whether it should be approved or requires revision.

## Input Data

### Proposed Annotation
**Cell Type**: <<cell_type>>

### Annotation Reasoning
<<draft_response_wreason>>

### Cluster Expression Profile
<<exp_summary_table>>

### Sample Metadata
<<metadata_table>>

### Reference Markers
<<cell_markers_table>>

## Evaluation Criteria

1. **Marker Consistency**: Do the proposed markers match the observed expression pattern?
2. **Biological Plausibility**: Is the cell type assignment biologically reasonable given the sample context?
3. **Specificity**: Is the annotation at the appropriate resolution level (not too broad or too narrow)?
4. **Evidence Quality**: Are there sufficient supporting markers, or is the annotation based on weak evidence?
5. **Contradiction Handling**: Has the annotator adequately addressed any conflicting evidence?

## Output Requirements

Provide your evaluation with detailed reasoning, addressing each criterion above. Then append your final decision in the following JSON format:

```json
{
    "final decision": "Approved"  // or "Disapproved"
}
```

If "Disapproved", explain what needs to be corrected and suggest a better annotation.

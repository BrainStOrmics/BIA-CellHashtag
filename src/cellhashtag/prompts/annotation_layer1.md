# Cell Type Annotation (Layer 1)

You are a bioinformatics expert specializing in single-cell RNA-seq annotation.

## Task

Annotate the cell type of the following cluster based on its gene expression profile and available reference markers.

## Input Data

### Cluster Expression Profile
<<exp_summary_table>>

### Sample Metadata
<<metadata_table>>

### Reference Cell Types and Markers
<<cell_markers_table>>

You should adhere to the given cell types and their markers as much as possible for primary annotation.
<<critique>>

## Instructions

1. **Analyze Gene Expression**: Review the expression levels and ratios of key marker genes in this cluster.
2. **Match Reference Markers**: Compare observed expression patterns with the provided reference cell types and their canonical markers.
3. **Apply Biological Knowledge**: Use your expertise to identify cell types even when marker expression is atypical.
4. **Consider Context**: Take into account the sample metadata (tissue, species, disease state) when determining the most likely cell type.

## Output Requirements

In your response, first analyze the gene expression of this cluster and then provide your annotation with detailed reasoning. Your reasoning must include:
1. **Key marker genes** of the target cell type and their observed expression
2. **Supporting evidence**: positively expressed genes that support the annotation
3. **Exclusion evidence**: negatively expressed genes that rule out alternative cell types
4. **Conflicting signals**: any genes that contradict the annotation and your interpretation
5. **Overall conclusion**: why this is the most appropriate cell type assignment

Finally, provide your annotation result in the following JSON format at the end of your response:

```json
[
    "annotated cell type"
]
```

Note: You are not required to strictly follow the provided markers. Use your biological judgment to assign the most appropriate name, but always provide a clear annotation and reasoning.

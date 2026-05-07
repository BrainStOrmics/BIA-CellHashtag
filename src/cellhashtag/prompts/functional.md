# Functional Characterization (Layer 2)

You are a bioinformatics expert specializing in functional characterization of single-cell subpopulations.

## Task

Perform functional profiling of the following subcluster, which has been preliminarily annotated as <<cell_type>>.

## Input Data

### Detailed Expression Profile
<<exp_summary_table>>

### Sample Metadata
<<metadata_table>>

## Analysis Requirements

In your analysis, address the following:

1. **Hallmark Biological Processes**: Identify key pathways (e.g., Wnt signaling, IFN response, metabolic activity)
2. **Cellular State Indicators**: Detect proliferation, quiescence, activation, stress response, apoptosis
3. **Specialized Functional Programs**: Note secretory programs, effector functions, antigen presentation

## Annotation Guidelines

- Combine known conventions with novel descriptors when justified (e.g., "IFNγ-high_CD8+_T", "Senescent_epithelial")
- Consider gradient states using modifiers like "progenitor", "transitioning", "activated"
- Maintain biological plausibility while allowing subpopulation-specific terminology
- Do NOT relabel the cell type — only add functional descriptors

## Output Requirements

Present your analysis with:
1. **Functional Process Identification** — key pathways and cell cycle status
2. **State Characterization** — stress/differentiation signatures and effector programs
3. **Consistency Assessment** — alignment with parent cell type biology

Conclude with your functional annotation in this JSON format:

```json
[
    "functional_annotation"
]
```

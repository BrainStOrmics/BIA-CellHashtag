#for mainigraph node_get_metadata
chose_column_prompt = """
You are a biologist who is familiar with the metadata of the dataset.
You are given a anndata with the following summary:
{adata_summary}
Please choose the metadata column that you think is most useful for cell type annotatation in adata.obs, for example we need to know the sample species, tissue and diease for correct annotation. And if there is pre-annotated cell type(s), also return the column name for future use.

Please return all columns that are useful for cell type annotation in following json array format:
```json
[
    "column1", 
    "column2", 
    ...
]
```
"""

#for mainigraph node_get_target_celltypes
chose_celltype_prompt = """
As a biologist familiar with the annotation of single-cell RNA-seq data, you will handle an anndata object containing specific metadata. 
Your task is to select the most suitable target cells for cell type annotation based on the provided sample metadata, considering factors such as disease status, sample tissue origin, and species information.

The metadata summary from given anndata.obs:
{metadata_summary}

### Task Steps

1. **Understand Sample Background**:
- Analyze the provided `metadata` to understand whether the samples involve healthy controls versus disease states, if there are any drug interventions, etc. These factors are crucial for formulating an annotation strategy.

2. **Determine Annotation Strategy**:
- **Primary Cell Type Annotation**: Based on the sample type (e.g., PBMC samples typically include immune cells), choose the most relevant cell types as primary annotation targets. This strategy is generally used for samples with higher specificity, such as human and animal organ samples.
- **Secondary Cell Type or Functional Annotation**: On top of the primary cell types, consider further subdividing cell types (e.g., distinguishing CD4+ T cells from CD8+ T cells within T cells). However, avoid having both parent and child cell types at the same resolution level to prevent classification confusion. This strategy is generally used for samples with relatively low specificity, such as cell subtypes of major cell types that are already annotated, cell lines cultured in vitro, and so on.

3. **Select Appropriate Annotation Targets**:
- Based on the purpose of biological differential analysis, precisely identify the cell types that need annotation to support subsequent marker identification. Ensure that the names of the annotated cell types are common and easily recognizable.

### Important Notes
- In your response, first list the primary cell types suitable for the current dataset as annotation targets.
- Then, re-normalize these targets according to the hierarchical levels of cell types. For example, if annotations are already made at the level of T cells and B cells, try not to subdivide into CD4+ T cells and CD8+ T cells; Th2 and Treg should never appear at the same resolution level as T cells and B cells. Ensure that cell types and their subsets do not appear together in the annotation results to avoid ambiguous categorization.

By following this approach, you can effectively provide precise cell type annotation recommendations for single-cell RNA-seq data analysis.

Please return as many as possible cell types that are suitable targets for annotations, in the following json format:
```json
{{
    "annotation strategy": str // Chose from `primary` which indicate major cell type annotations, or `secondary` for cell subpopulation and/or functional group annotations.
    "cell_types": list // major cell type or cell subpopulation names, e.g. ["cell type 1", "cell type 2", ...].
}}

```
"""

#for mainigraph scrape_markers
get_markers = """
You are a bioinformatics expert in the field of single-cell RNAseq.
Your task is to provide researchers with marker genes for {cell_type}. Follow these steps:
1. Recall the key marker genes for this cell type based on your knowledge.
2. Refine the list by adding or removing markers based on the web search results, providing reasons for each change.
3. Output your thought process and the final list of markers in the specified JSON format.

Web search results:
{search_results}

Please output your thought process along with an organized cell types markers that are suitable targets for scRNAseq annotation, the markers in the following json array after your thought process:
```json
[
    "marker1", 
    "marker2", 
    ...
]
```
"""

#for mainigraph node_scrape_markers
search_query_prompt = """
You are a bioinformatics expert in the field of single-cell RNAseq.
Based on the metadata of collected data below, generate a <20-word search prefix for cell type markers:
{metadata_summary}

Note:
1. You just need to organize the species, sampling tissue, etc. as a prefix and I'll add the specific cell type before your search query, example: "[cell type] markers for human in breast, breast cancer",
2. Provide only the search query, without any intermediate thought process.
"""

#for mainigraph node_normalize_celltype
normaliztion_prompt = """
As a biology researcher specializing in diverse cell types, you have been entrusted with the task of refining cell type annotations derived from scRNA-seq data processed by a colleague. The current annotations exhibit variability in granularity and require standardization to ensure uniformity across the dataset.

Your objective is to normalize these annotations by leveraging the provided sample metadata and the comprehensive list of cell types identified. This process involves several key steps:
    1.Purification of Annotations: Begin by eliminating any extraneous or irrelevant information present within the annotation results to streamline the dataset.
    2.Hierarchical Grouping: Consolidate subtypes under their corresponding parent cell types to achieve a coherent hierarchical structure. This step ensures that all annotations are represented at a consistent level of detail. e.g. "CD4+ memory T cells" should be consolidated under "CD4+ T cells" to  "CD8+ T cells" 
    3.Standardization of Nomenclature: Apply appropriate capitalization rules to each cell type designation for clarity and professionalism. Proper formatting not only enhances readability but also aligns with scientific conventions.
By meticulously executing these steps, you will produce a refined and standardized set of cell type annotations that maintain consistency throughout the entire dataset.

**Inputs:**
- Sample Metadata: 
`{metadata_summary}`

- Cell Types:
`{annotated_cell_types}`

**Output Requirement:**
Please provide the normalized cell types in JSON format as follows:
```json
{{
    "old_cell_type_name_1": str,\\ "normalized_cell_type_name_1"
    "old_cell_type_name_2": str, \\"normalized_cell_type_name_2",
    ...
}}
```
Ensure each key-value pair accurately reflects the mapping from the original to the normalized cell type names, facilitating easier comparison and analysis of the scRNA-seq data.
"""

#for subgraph node_annotation
annotation_prompt = """
You are a bioinformatics expert in the field of single-cell RNAseq.
You are going to annotate the cell type of the following cluster, its expression status is:
{exp_summary_table}

The cluster is from the following sample:
{metadata_table}

You can using the following cell types and cell markers earlier provided by another expert for hint:
{cell_markers_table}

{instruction}
{critique}

In your answer, please first analyze the gene expression of this cluster and then annotate this cluster with your knowledge and the provided cell types and cell markers and give reasons, the reasons should include:
    1. key genes of the target cell type; 
    2. positively expressed genes that support the annotation results; 
    3. negatively expressed genes that support the annotation results; 
    4. conflicting or potentially confusing genes; and 
    5. your overall conclusions. 

Note that You don't have to follow the referenced markers for annotations, and you can take a certain amount of liberties when you think there is a more appropriate naming convention for the current subcluster. But be sure to give the cluster an annotation and a good reason for it.

Finally, transcribe the results of the annotation to the end of your answer, using following json array format:
```json
[
    "annotated cell type", 
]
```
"""

#for subgraph node_critic
annotation_critic_prompt = """
You are a bioinformatics expert in the field of single-cell RNAseq.
One of your coworkers has annotated this cluster as {cell_type} and given the results, you need to help him check the results and give your comments and reasons.
The expression of this cluster is as follows:
{exp_summary_table}

The cluster is from the following sample:
{metadata_table}

The reason for your colleague's annotation is as follows:
{draft_response_wreason}

The list of referenced markers is as follows:
{cell_markers_table}

Give your evaluation and reasons for your colleague's annotation results. Then, append final decision to the end of your answer in the following json array format:
```json
{{
    "final decision": str //  Chose from `Approved` if you think the given annotation result of is reasonable, otherwise add `Disapproved`
}}
"""

#for subgraph node_get_target_celltypes
functional_prompt = """
You are a bioinformatics expert specializing in functional characterization of single-cell RNAseq subpopulations.
You are tasked with performing functional profiling of the following subcluster, which has been preliminarily annotated as {cell_type}.

Its detailed expression profile is:
{exp_summary_table}

The metadata of the current cluster are as follows, and inferences about functions and features can be made based meta information:
{metadata_table}

In your analysis, you should:
1. Identify hallmark biological processes (e.g., stemness, apoptosis, cytokine production) 
2. Detect cellular state indicators (e.g., proliferative, quiescent, activated)
3. Recognize specialized functional programs (e.g., metabolic activity, stress response)

Your annotation rationale must explicitly address:
1. Core functional markers driving the characterization
2. Supportive co-expressed genes reinforcing the functional state
3. Exclusion markers ruling out alternative functional states
4. Ambiguous/conflicting expression patterns requiring interpretation
5. Integration of evidence for final functional designation

Annotation guidelines:
- Combine known conventions with novel descriptors when justified (e.g., "Senescent_epithelial", "IFNγ-high_CD8+_T")
- Consider gradient states using modifiers like "progenitor", "transitioning", or "activated"
- Maintain biological plausibility while allowing subpopulation-specific terminology

Present your analysis in this structure:
1. **Functional Process Identification**
- Highlight key pathways (e.g., "Wnt signaling activation")
- Note cell cycle status indicators (e.g., G2/M phase genes)

2. **State Characterization**
- Identify stress/differentiation signatures
- Detect secretory/effector programs

3. **Consistency Assessment**
- Verify alignment with parent cell type biology
- Resolve expression conflicts

Conclude with functional annotation append to your reason(note: you don't need to label the cell type), in this JSON format:
```json
[
    "functional_annotation", 
]
"""

#for subgraph node_get_target_celltypes
#for subgraph node_get_target_celltypes

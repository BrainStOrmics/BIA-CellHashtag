"""
Utility functions and classes.

"""
#import types
from pandas import DataFrame 
from anndata import AnnData
from langgraph.graph.state import CompiledStateGraph
#import packages
import json 
import scanpy as sc
import numpy as np
import pandas as pd
from io import StringIO
from contextlib import redirect_stdout

#draw graph
from IPython.display import Image, display

def call_chatllm_openai(api_key, api_base, model_name):
    llm = ChatOpenAI(
        openai_api_key = api_key,
        openai_api_base=api_base,
        model = model_name)
    return llm

def parse_jsonfromcontent(content: str):
    '''
    Extracts and parses JSON data from a string containing JSON data.
    This function looks for '```json'' and '```'' tags in the input string and extracts the JSON data between them.
    It then tries to parse the extracted string into Python objects (usually a dictionary or list).

    Parameters: 
    content (str): A string containing the JSON data.

    Returns: 
    dict or list: the parsed Python object.

    Raises: 
    ValueError: If the '```json'' or '```' token is not found, or JSON decoding fails
    '''
    start_marker = '```json'
    end_marker = '```'
    # search JSON start position
    start = content.find(start_marker)
    if start == -1:
        raise ValueError("Cannot find JSON start string: '```json'")
    # calculate JSON data start position
    start += len(start_marker)
    # search JSON end position
    end = content.find(end_marker, start)
    if end == -1:
        raise ValueError("Cannot find JSON end string '```'")
    # extract JSON string and remove whitespace     
    json_str = content[start:end].strip()
    # trying to parse JSON data
    try:
        return json.loads(json_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Failed to parse JSON: {e}") 

def get_cluster_hvgs(
        adata: AnnData, 
        groupby='leiden', 
        n_genes=5, 
        plot=False) -> list | None:
    """
    Get highly variable genes (HVGs) for each group in the AnnData object.
    
    This function checks if 'rank_genes_groups' is already computed in the AnnData object.
    If not, it computes it using the specified `groupby` column. It then extracts the top
    `n_genes` HVGs for each group and optionally plots a dotplot if requested.
    
    Parameters:
        adata (AnnData): The AnnData object containing gene expression data.
        groupby (str): Column name in `adata.obs` to group by. Default is 'leiden'.
        n_genes (int): Number of top genes to extract per group. Default is 5.
        plot (bool): Whether to generate a dotplot of ranked genes. Default is False.
    
    Returns:
        list: A list of HVGs from all groups.
    
    Raises:
        ValueError: If `groupby` is not a column in `adata.obs`.
    """
    # Validate that the groupby column exists in adata.obs
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs")
    # Check if ranked genes are pre-computed; if not, compute them
    if 'rank_genes_groups' not in adata.uns.keys():
        sc.tl.rank_genes_groups(adata, groupby=groupby)
    # Optionally plot a dotplot of ranked genes
    if plot:
        sc.pl.rank_genes_groups_dotplot(
            adata, 
            groupby=groupby, 
            n_genes=25,
            save='HVGfor'+groupby+'.pdf')
    # Get unique group labels (e.g., clusters)
    groups = adata.obs[groupby].unique().tolist()
    # Collect top n_genes HVGs for each group
    hvgs = []
    for group in groups:
        hvgs.extend(adata.uns['rank_genes_groups']['names'][group][:n_genes])
    return hvgs

def df2markdownTable(
        df:DataFrame) -> str| None:
    """
    Convert a pandas DataFrame to a markdown table format.

    This function takes a pandas DataFrame and returns a string representing the data
    in a properly formatted markdown table, including a header row, a separator row,
    and data rows.

    Parameters:
        df (pd.DataFrame): The input DataFrame to convert into a markdown table.

    Returns:
        str: A string containing the markdown table representation of the DataFrame.

    Raises:
        ValueError: If the input is not a pandas DataFrame.
    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame")
    # Create the header row with column names separated by ' | '
    header = '| ' + ' | '.join(df.columns) + ' |'
    # Create the separator row with '---' for each column, required for markdown formatting
    separator = '| ' + ' | '.join(['---'] * len(df.columns)) + ' |'
    # Build the data rows by iterating over the DataFrame
    data_rows = []
    for _, row in df.iterrows():
        # Convert each value in the row to a string and join with ' | '
        row_str = '| ' + ' | '.join(str(value) for value in row) + ' |'
        data_rows.append(row_str)
    # Combine header, separator, and data rows into a single string with newlines
    table = header + '\n' + separator + '\n' + '\n'.join(data_rows)
    return table

def get_exp_summary(
        adata:AnnData, 
        genes:list, 
        use_raw =False, 
        exp_cutoff = 0.01) -> str | None:
    """
    Generate a Markdown-formatted table summarizing gene expression levels and ratios from an AnnData object.
    
    The expression level is computed as the average expression value of a gene across all cells.
    The expression ratio is the percentage of cells where the expression exceeds a given threshold.
    
    Parameters:
        adata (AnnData): The AnnData object with gene expression data.
        genes (list): List of gene names to summarize.
        use_raw (bool): If True, use raw data (adata.raw); otherwise, use processed data. Default is False.
        exp_cutoff (float): Threshold for expression ratio calculation. Default is 0.01.
    Returns:
        str: A Markdown-formatted table with expression summaries.
    Raises:
        ValueError: If a gene is not uniquely identified in the AnnData object.
    """
    if use_raw and adata.raw is not None:
        adata_exp = adata.raw.to_adata()
    else:
        adata_exp = adata.copy()
    #calculate expression level and ratio
    expression_summary = '==========\n'
    expression_summary += '| gene name | expression level | expression ratio (%) |\n'
    expression_summary += '|-----------|------------------|----------------------|\n'
    # Analyze each gene
    for gene in genes:
        if gene not in adata_exp.var_names:
            continue #skip this gene
        else:
            exp_matrix = adata_exp[:,adata_exp.var_names == gene].X
            if exp_matrix.shape[1] != 1:
                raise ValueError(f"Gene '{gene}' is not unique in the AnnData object.")
            if hasattr(exp_matrix, 'toarray'):
                exp_matrix = exp_matrix.toarray()
            # Compute expression metrics
            try:
                exp_level = np.mean(exp_matrix)
                exp_ratio = adata_exp[exp_matrix>exp_cutoff].n_obs/adata_exp.n_obs * 100
            except:
                exp_level = 0
                exp_ratio = 0
            #convert to markdown
            expression_summary+= f'| {gene} | {exp_level:.2f} | {exp_ratio:.2f} |\n'
    # Complete the table
    expression_summary+= '==========\n'
    return expression_summary

def data_perception(data):
    """
    Analyzes basic properties and a summary of the given data.
    
    Parameters:
    - data: The input data object to be analyzed.
    
    Returns:
    - A tuple containing the type of the input data and its string representation.
    """
    #get datatype
    data_type = type(data)
    #get summary 
    capture = StringIO()
    with redirect_stdout(capture):
        print(data)
    data_summary = capture.getvalue()
    capture.close()
    return data_type, data_summary

def get_dfcol_summaries(
        df: pd.DataFrame, 
        columns: list) -> str:
    """
    Generates summaries of value counts for specified columns in a pandas DataFrame.
    
    Parameters:
    - df: The pandas DataFrame containing the data.
    - columns: A list of column names for which to generate value count summaries.
    
    Returns:
    - A formatted string containing the value counts for each specified column.
    """
    summaries = '============\n'
    for col in columns:
        if col not in df.columns:
            summaries += f"Column '{col}' does not exist in the DataFrame.\n\n"
            continue
        # Get value counts for the column
        col_summary = df[col].value_counts().to_string()
        # Format the summary
        summaries += f'Value counts of column {col}:\n{col_summary}\n\n'
    summaries += '============\n'
    return summaries

def set_web_scraper(
        api_key, 
        scraper = 'Tavily',
        k_web_per_type=7,
        ):
    if scraper == 'Tavily':
        from langchain_community.retrievers import TavilySearchAPIRetriever
        web_retriever = TavilySearchAPIRetriever(k=k_web_per_type, api_key = api_key)
    return web_retriever

def draw_graph(graph:CompiledStateGraph):
    display(Image(graph.get_graph(xray=1).draw_mermaid_png()))

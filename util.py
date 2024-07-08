"""
Utility functions for GRN analysis and visualization.

This module contains a collection of functions used for preprocessing RNA-seq data,
inferring gene regulatory networks, performing Gene Set Enrichment Analysis (GSEA),
and analyzing isoform expression patterns across different tissues.

1. Data Preprocessing:
   - extract_and_save_gtex_columns: Extract GTEX columns from input files to make a sample ID list
   - read_transcript, read_gene: Read and process transcript and gene data for each tissue
   - protein_filtering, protein_gene_filtering: Filter for protein-coding genes
   - filtering_low_expression: Remove low-expression genes

2. GRN Inference:
   - run_grnboost: Run GRNBoost2 algorithm for GRN inference
   - grn_mapping: Construct an alternative splicing-aware GRN

3. Differential Expression Analysis:
   - gene_expression_sum: Combine gene expression data from each tissue
   - meta_data: Create metadata for differential expression analysis

4. Isoform Analysis:
   - isoform: Calculate isoform expression patterns between tissues
   - merge_grn_with_percentages: Merge GRN data with isoform percentages
   
5. GSEA:
   - map_names_GSEA_target: Prepare data for GSEA
   - run_gsea: Perform Gene Set Enrichment Analysis


Usage:
This module is intended to be imported and used in Jupyter notebooks or other
Python scripts for GRN analysis.

Author: Hyunjung Wang
Date: 7/7/2024

"""

import gzip
import pandas as pd
import os
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
import pandas as pd
from pybiomart import Server
import gseapy as gp
from gseapy import dotplot
import matplotlib.pyplot as plt
import sys
import numpy as np

from distributed import Client, LocalCluster
import dask.dataframe as dd

# Define the relative path to the parent directory
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
# Add the parent directory to sys.path
sys.path.append(parent_dir)

tissue_gene_paths = {
    "bladder": "bladder_tr_gene.txt",
    "uterus": "uterus_tr_gene.txt",
    "kidney": "kidney_tr_gene.txt",
    "brain": "brain_tr_gene.txt",
    "cervix_ectocervix": "cervix_ectocervix_tr_gene.txt",
    "fallopian_tube":"fallopian_tube_tr_gene.txt"
}

tissue_sample_paths = {
    "bladder": "bladder_samples.txt",
    "uterus": "uterus_samples.txt",
    "kidney": "kidney_samples.txt",
    "brain": "brain_samples.txt",
    "cervix_ectocervix": "cervix_ectocervix_samples.txt",
    "fallopian_tube":"fallopian_tube_samples.txt"
}

def extract_and_save_gtex_columns(input_file, output_file):
    """
    Extract columns starting with 'GTEX' from a given file and save them to a new file.
    This function creates a GTEX sampleIDs list.
    
    Parameters:
    input_file (str): The path to the input file.
    output_file (str): The path to the output file.
    """
    # Read the file while skipping the first two rows
    df = pd.read_csv(input_file, sep='\t', skiprows=2)
    
    # Extract columns that start with 'GTEX'
    gtex_columns = [col for col in df.columns if col.startswith('GTEX')]
        
    # Create a DataFrame with the GTEX column names
    gtex_columns_df = pd.DataFrame(gtex_columns, columns=["GTEX_Column_Names"])
    
    # Save the GTEX column names to a new CSV file without header
    gtex_columns_df.to_csv(output_file, index=False, header=False)

    
def read_transcript(gz_file_path, sample_path, output_file_path):
    """
    Read transcript data from a gzipped file, select specific columns, and save to an output file.

    Args:
        gz_file_path (str): Path to the gzipped input file.
        sample_path (str): Path to the file containing selected column names.
        output_file_path (str): Path to the output file where processed data will be saved.
    """
    # Load the selected columns(samples)
    selected_columns = pd.read_csv(sample_path, header=None).squeeze()

    # Open the GZ file and save 'transcript_id' and selected columns to the output file
    with gzip.open(gz_file_path, 'rt') as gz_file, open(output_file_path, 'w') as output_file:
        # Skip the first two rows
        next(gz_file)
        next(gz_file)
        
        # Read the header to get column indices
        header = gz_file.readline().strip().split('\t')
        column_indices = [header.index(col) for col in selected_columns]

        # Write the header for 'transcript_id' and selected columns to the output file
        output_file.write("gene_id\t"+ "transcript_id\t" + '\t'.join(selected_columns) + '\n')

        # Write 'transcript_id' and selected columns to the output file for the remaining rows
        for line in gz_file:
            # Split the line into columns based on the tab separator
            columns = line.strip().split('\t')

            # Extract 'transcript_id' and selected columns using column indices
            gene_id = columns[1]
            transcript_id = columns[0]

            selected_values = [columns[idx] for idx in column_indices]

            # Write 'transcript_id' and selected columns to the output file
            output_file.write("{}\t{}\t{}\n".format(gene_id, transcript_id, '\t'.join(selected_values)))

def read_gene(gz_file_path, sample_path, output_file_path):
    """
    Read gene data from a gzipped file, select specific columns, and save to an output file.

    Args:
        gz_file_path (str): Path to the gzipped input file.
        sample_path (str): Path to the file containing selected column names.
        output_file_path (str): Path to the output file where processed data will be saved.
    """
    # Load the selected columns(samples)
    selected_columns = pd.read_csv(sample_path, header=None).squeeze()

    # Open the GZ file and save 'transcript_id' and selected columns to the output file
    with gzip.open(gz_file_path, 'rt') as gz_file, open(output_file_path, 'w') as output_file:
        # Skip the first two rows
        next(gz_file)
        next(gz_file)
        
        # Read the header to get column indices
        header = gz_file.readline().strip().split('\t')
        column_indices = [header.index(col) for col in selected_columns]

        # Write the header for 'transcript_id' and selected columns to the output file
        output_file.write("gene_id\t"+ '\t'.join(selected_columns) + '\n')

        # Write 'transcript_id' and selected columns to the output file for the remaining rows
        for line in gz_file:
            # Split the line into columns based on the tab separator
            columns = line.strip().split('\t')

            # Extract 'transcript_id' and selected columns using column indices
           # Description = columns[1]
            Name = columns[0]

            selected_values = [columns[idx] for idx in column_indices]

            # Write 'transcript_id' and selected columns to the output file
            output_file.write("{}\t{}\n".format(Name, '\t'.join(selected_values)))



def protein_filtering(df):
    """
    Filter the input DataFrame to keep only protein-coding transcripts.

    Args:
        df (pd.DataFrame): Input DataFrame containing transcript expression data.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only protein-coding transcripts.
    """
    # filter lastest version of stable id
    def get_latest_stable_id(df):
        transcript_id, version = df.iloc[0, 1].split('.')
        return df[df['transcript_id'].str.split('.').str[0] == transcript_id]  # Filter by gene ID

    # Read protein_coding.txt into df_protein DataFrame
    protein_coding_path = os.path.join(parent_dir, "protein_coding.txt")

    df_protein = pd.read_csv(protein_coding_path)
    df = df.groupby('transcript_id').apply(get_latest_stable_id).reset_index(drop=True)

    df['gene_id'] = df['gene_id'].str.split('.').str[0]
    df['transcript_id'] = df['transcript_id'].str.split('.').str[0]
    
    # Merge based on 'gene_id' and 'transcript_id', drop 'Gene type' column
    merged_df = pd.merge(df, df_protein.drop(columns=['Gene type']), 
                         left_on=['gene_id', 'transcript_id'], 
                         right_on=['Gene stable ID', 'Transcript stable ID'], 
                         how='inner')
    
    # Drop 'Gene stable ID version' and 'Transcript stable ID version' columns
    merged_df = merged_df.drop(columns=['Gene stable ID', 'Transcript stable ID'])

    merged_df = merged_df.drop(columns=['gene_id'])
    
    # Transpose the DataFrame without keeping the index
    df_transposed = merged_df.transpose().reset_index(drop=True)

    # Set the header to be the gene name
    df_transposed.columns = df_transposed.iloc[0]
    df_transposed = df_transposed.drop(0)

    
    return df_transposed


def protein_gene_filtering(tissue, gene_data):
    """
    Filter gene data based on protein coding information and merges it with transcript expression data.

    Args:
        filtered_genes_file_path (str): Path to the filtered gene data file.
        gene_data (str): Path to the gene data.

    Returns:
        pandas.DataFrame: The merged DataFrame containing filtered genes and tissue data.
    """
    # Define the relative path to the parent directory
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
    # Add the parent directory to sys.path
    sys.path.append(parent_dir)

    def get_latest_stable_id(df):
        gene_id, version = df.iloc[0, 0].split('.')
        return df[df['gene_id'].str.split('.').str[0] == gene_id]  # Filter by gene ID

    # Read the filtered gene data

    df = pd.read_csv(gene_data, sep='\t')
    df = df.groupby('gene_id').apply(get_latest_stable_id).reset_index(drop=True)
    df['gene_id'] = df['gene_id'].str.split('.').str[0]
  
    # Construct the relative path to protein_coding_genes.txt
    protein_genes_path = os.path.join(parent_dir, "protein_coding_genes.txt")

    df_protein = pd.read_csv(protein_genes_path)
    
    # Merge based on 'gene_id' and 'transcript_id', drop 'Gene type' column
    merged_df = pd.merge(df, df_protein.drop(columns=['Gene type']), 
                         left_on=['gene_id'], 
                         right_on=['Gene stable ID'], 
                         how='inner')
    
    # Drop 'Gene stable ID version' and 'Transcript stable ID version' columns
    merged_df = merged_df.drop(columns=['Gene stable ID','Transcript stable ID'])
    
    
    df_transposed = merged_df.transpose().reset_index(drop=True)
    headers = df_transposed.iloc[0].values
    df_transposed.columns = headers
    df_transposed.drop(index=0, axis=0, inplace=True)
    df_transposed = df_transposed.reset_index(drop=True)

    # Read transcript data
    tissue_df = pd.read_csv(f"{tissue}_filtered_tr.txt", index_col=0)

    # Concatenate tissue data with filtered gene data
    merge_df2 = pd.concat([tissue_df, df_transposed], axis=1)

    # Save merged DataFrame to file
  #  merge_df2.to_csv(f"{tissue}_gene_filtered.txt")
    
    return merge_df2



def filtering_low_expression(df, threshold_fraction, min_value):
    """
    Filter out low-expression genes in RNA-seq data.

    Parameters:
    - df: DataFrame containing gene expression data
    - threshold_fraction: Fraction of the total sequence length used as the threshold for filtering
    - min_value: Minimum value required for a gene to be considered expressed

    Returns:
    - df_filtered: Filtered DataFrame with low-expression genes removed
    """
    # Calculate the threshold based on the total sequence length
    threshold = len(df) * threshold_fraction

    # Remove genes where more than the threshold fraction of the values are zero
    zero_counts = df.eq(0).sum()
    genes_to_remove = zero_counts[zero_counts > threshold].index
    df_filtered = df.drop(columns=genes_to_remove)

    # Filter the DataFrame to keep only genes where all values are greater than the minimum value
    genes_to_keep = (df_filtered > min_value).all()
    df_filtered = df_filtered.loc[:, genes_to_keep]

    return df_filtered
    



def preprocessing_data(tissue, available_tissues, threshold_fraction, min_value):
    """
    Preprocess raw data for a specific tissue.
    This function performs protein-coding gene filtering, removal of low expression genes, and data formatting.
    It's a fundamental step in preparing the raw data for downstream GRN inference.

    Args:
        tissue (str): Name of the tissue to preprocess.
        available_tissues (list): List of available tissue names for validation.

    Returns:
        None. Saves preprocessed data to various output files.
    """
    
    if tissue not in available_tissues not in available_tissues:
        print("Invalid tissues selected. Please choose from available_tissues.")
        return

    sample_tissue = f"{tissue}_samples.txt"

    transcript_tissue = f"transcript_{tissue}.txt"

    # Read and process tissue transcript data
    df_tissue = pd.read_csv(transcript_tissue,sep='\t')
    #protein-coding filtering  
    df_tissue_protein= protein_filtering(df_tissue)
    df_tissue_protein.to_csv(f"{tissue}_protein_filtered.txt",index=False)
    df_tissue_protein=pd.read_csv(f"{tissue}_protein_filtered.txt")

    #low expression filtering
    df_tissue_filtered = filtering_low_expression(df_tissue_protein, threshold_fraction, min_value)
    df_tissue_filtered.to_csv(f"{tissue}_filtered_tr.txt")

    # Get unique column values for tissue-.
    df_tissue_unique = df_tissue_filtered.columns.drop_duplicates()
    df_tissue_unique.to_series().to_csv(f"{tissue}_column.txt", index=False, header=False, sep='\t')
    

    # Read and process tissue gene data
    gene_tissue = f"gene_{tissue}.txt"
    #protein-coding filtering and merging with transcript data
    tissue_protein_filtered = protein_gene_filtering(tissue, gene_tissue)

    # low expression filtering
    tissue_gene_filtered = filtering_low_expression(tissue_protein_filtered,threshold_fraction, min_value)
    # Save the filtered data
    tissue_gene_filtered.to_csv(f"{tissue}_tr_gene.txt")


def tf_names(columns,output_tf_id):
    """
    Filter transcription factor IDs based on a given list and save the results.

    Args:
        columns (str): Path to the file containing the list of transcription factor IDs.
        output_tf_id (str): Path to save the filtered transcription factor IDs.
    """
    # Define the relative path to the parent directory
    transcription_factor = os.path.join(parent_dir, "TF_list.txt")

   # Add the parent directory to sys.path
    sys.path.append(parent_dir)
    with open(columns, "r") as tf_file:
        tf_id = set(line.strip() for line in tf_file)

    # Initialize a list to store filtered data
    filtered_data = []

    # Read and filter data from "TF_list.txt"
    with open(transcription_factor, "r") as annot_file:
        # Assuming the first line contains the header
        header = annot_file.readline().strip().split(",")

        # Find the indices of relevant columns
        ensembl_transcript_id_index = header.index("ensembl_transcript_id")
        external_gene_name_index = header.index("external_gene_name")

        # Add the header to filtered_data
        filtered_data.append(header)

        # Iterate through the remaining lines in the file
        for line in annot_file:
            # Split the line into values based on comma separator
            values = line.strip().split(",")

            # Extract transcript ID and gene name from the line
            ensembl_transcript_id = values[ensembl_transcript_id_index]

            # Check if the transcript ID is in the set of TF IDs
            if ensembl_transcript_id in tf_id:
                # If yes, add the line to the filtered data
                filtered_data.append(values)

    # Write the filtered data to a new file named "kidney_tf_id.txt"
    with open(output_tf_id, "w") as filtered_file:
        # Write each line of filtered data to the file
        for line in filtered_data:
            filtered_file.write(line[0] + "\n")
            #filtered_file.write(",".join(line) + "\n")


def grnboost(df,tf,output):
    """
    Run the GRNBoost2 algorithm to infer gene regulatory networks.

    Args:
        df (pd.DataFrame): Expression data DataFrame.
        tf (str): Path to the file containing transcription factor names.
        output (str): Path to save the output network.

    Returns:
        None. Saves the inferred network to the specified output file.
    """

    # Reads each line of the file, with each line containing the name of a TF.
    tf_names = load_tf_names(tf)

    # Run GRNBoost2
    network = grnboost2(expression_data=df,tf_names=tf_names)

    # Save the network
    network.to_csv(output, sep='\t', header=True, index=False)



def run_grnboost(tissue):
    """
    Run the GRNBoost2 algorithm for a specific tissue.

    Args:
        tissue (str): Name of the tissue to analyze.

    Returns:
        None. Processes the data and saves the results.
    """
    # Process tf_names 
    tf_names(f"{tissue}_column.txt", f"{tissue}_tf_id.txt")

    # Load the filtered data 
    tissue_filtered = pd.read_csv(f"{tissue}_tr_gene.txt",index_col=0)
   
    # Set up a local cluster for parallel processing
    cluster = LocalCluster()
    # Create a client to interact with the cluster

    client = Client(cluster)
    # Run GRNBoost
    grnboost(tissue_filtered, f"{tissue}_tf_id.txt", f"{tissue}_network.tsv")

    # Shut down Dask client
    client.close()
    cluster.close()

def grn_mapping(tissue):
    """
    Constructing an alternative splicing-aware GRN structure using Ensemble BioMart.
    If the target is a transcript, map transcript ID to gene ID
    else map to 0
    Map TF transcript IDs to gene IDs

    Args:
        tissue (str): Name of the tissue.

    Returns:
        pd.DataFrame: Mapped and sorted network data.
    """
    # Read the network data
    network_data = pd.read_csv(f'{tissue}_network.tsv', sep='\t')

    # Sort the DataFrame by 'importance' column in descending order
    network_data_sorted = network_data.sort_values(by='importance', ascending=False)

    # Connect to Ensembl Biomart server
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_gene_id', 'ensembl_transcript_id']

    # Query Biomart database to retrieve transcript ID to gene ID mapping
    result_mapping = dataset.query(attributes=attributes)

    # Convert result mapping to a dictionary for easy lookup
    transcript_to_gene_mapping = dict(zip(result_mapping['Transcript stable ID'], result_mapping['Gene stable ID']))

    # Define a function to get gene ID based on transcript ID
    def get_gene_id(target):
        if 'ENST' in target:  # If target contains 'ENST', it is a transcript ID
            return transcript_to_gene_mapping.get(target, 0)  # Get gene ID if found, otherwise return None
        else:
            return 0  # Return None for non-transcript IDs

    # Add a new column 'target_transcript' containing gene IDs corresponding to 'target' transcript IDs
    network_data_sorted['target_transcript_id'] = network_data_sorted['target'].apply(get_gene_id)

    # Iterate over the DataFrame rows
    for index, row in network_data_sorted.iterrows():
        # Check if the 'target' value starts with 'ENST'
        if row['target'].startswith('ENST'):
            # Exchange the values of 'target' and 'target_transcript'
            network_data_sorted.at[index, 'target'], network_data_sorted.at[index, 'target_transcript_id'] = row['target_transcript_id'], row['target']

    combined_df = pd.merge(network_data_sorted, result_mapping, left_on='TF', right_on='Transcript stable ID', how='inner')

    
    # Rename columns from the 'result' dataframe
    combined_df.rename(columns={'Gene stable ID': 'TF_gene', 'Transcript stable ID': 'TF_transcript_id','target':'target_gene'}, inplace=True)

    combined_df['tissue'] = tissue

    # Drop unnecessary columns
    combined_df = combined_df.drop(['TF'], axis=1)

    column_order = [ 'tissue', 'TF_gene', 'TF_transcript_id', 'target_gene','target_transcript_id', 'importance']

    combined_df = combined_df[column_order]

    # Save the modified dataframe to a CSV file
    combined_df.to_csv(f"{tissue}_df.csv", index=False)

    return combined_df


def isoform_gene_expression_data(data, samples,tissue):
    """
    Prepare expression data for a specific tissue to calculate tissue-specific isoform percentage within the same gene.

    Args:
        data (str): Path to the expression data file.
        samples (str): Path to the file containing sample names.
        tissue (str): Name of the tissue.

    Returns:
        pd.DataFrame: Processed expression data.
    """
    # Read the sample names from the samples file
    sample_names = pd.read_csv(samples, header=None, names=['index_column'])

    # Read the data file
    df = pd.read_csv(data)  # Assuming tab-separated data, adjust delimiter as needed
    df = df.reset_index(drop=True)
    df = df.drop(df.columns[0], axis=1)
    df.set_index(sample_names['index_column'], inplace=True)

    # Additional processing to merge transcript IDs and gene IDs
    transcript_ids = df.columns.tolist()

    # Connect to the Ensembl database
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_gene_id', 'ensembl_transcript_id']

    result = dataset.query(attributes=attributes)
    
    # Filter the result to keep only the rows corresponding to transcript IDs in the DataFrame
    filtered_result = result[result['Transcript stable ID'].isin(transcript_ids)]
    filtered_result.set_index('Transcript stable ID', inplace=True)

    # Reorder based on transcript IDs
    filtered_result = filtered_result.T

    # Reorder the DataFrame based on the index of filtered_result
    df_reordered = df.reindex(columns=filtered_result.columns)

    # Concatenate the DataFrames
    merged_df = pd.concat([filtered_result, df_reordered])
        # Remove the row with index name 'Transcript stable ID'
    merged_df = merged_df.drop('Transcript stable ID', axis=0, errors='ignore')

   
    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(f"{tissue}_expression.csv")

    return merged_df



 
def gene_expression_sum(tissue1, tissue2):
    """
    Combine gene expression data for two tissues and sum them up if they have duplicated genes.

    Args:
        tissue1 (str): Name of the first tissue.
        tissue2 (str): Name of the second tissue.

    Returns:
        pd.DataFrame: Combined gene expression data for both tissues.
    """
   
    # Load sample names
    tissue1_samples = pd.read_csv(tissue_sample_paths[tissue1],  header=None, names=['Gene stable ID'])
    tissue2_samples = pd.read_csv(tissue_sample_paths[tissue2],  header=None, names=['Gene stable ID'])
    
    # Load tissue data
    tissue1_df = pd.read_csv(tissue_gene_paths[tissue1])
    tissue2_df = pd.read_csv(tissue_gene_paths[tissue2])

    tissue1_df.reset_index(drop=True)
    tissue2_df.reset_index(drop=True)
    
    tissue1_df = tissue1_df.drop(tissue1_df.columns[0], axis=1)
    tissue2_df = tissue2_df.drop(tissue2_df.columns[0], axis=1)



    # Set index to Gene stable ID from sample names
    tissue1_df.set_index(tissue1_samples['Gene stable ID'], inplace=True)
    tissue2_df.set_index(tissue2_samples['Gene stable ID'], inplace=True)


    
    # Transpose the dataframes
    tissue1_df_T = tissue1_df.T
    tissue2_df_T = tissue2_df.T

    # Group by gene ID and sum expression values
    grouped_tissue1 = tissue1_df_T.groupby(tissue1_df_T.index.str.split('.').str[0]).sum()
    grouped_tissue2 = tissue2_df_T.groupby(tissue2_df_T.index.str.split('.').str[0]).sum()

    # Rename the index for each dataframe
    grouped_tissue1.rename_axis('Gene_ID', inplace=True)
    grouped_tissue1.columns.name = None
    
    grouped_tissue2.rename_axis('Gene_ID', inplace=True)
    grouped_tissue2.columns.name = None

    # Merge the dataframes based on gene ID
    combined_df = pd.merge(grouped_tissue1, grouped_tissue2, on='Gene_ID', how='outer')

    # Fill missing values with 0
    combined_df.fillna(0, inplace=True)

    # Save the combined dataframe to a file
    combined_df.to_csv(f"../differential_expression_analysis/{tissue1}_{tissue2}_gene_expression_sum.csv")


    return combined_df  
    
def meta_data(tissue1, tissue2):
    """
    Create a metadata file that maps sample IDs to their corresponding tissues for two tissues.
    Args:
        tissue1 (str): Name of the first tissue.
        tissue2 (str): Name of the second tissue.

    Returns:
        None. Saves the metadata to a CSV file.
    """
    # Load sample names without headers
    samples_tissue1 = pd.read_csv(f'{tissue1}_samples.txt', header=None, names=['Sample_ID'])
    samples_tissue2 = pd.read_csv(f'{tissue2}_samples.txt', header=None, names=['Sample_ID'])

    # Concatenate the sample DataFrames
    df = pd.concat([samples_tissue1, samples_tissue2], ignore_index=True)

    # Add a 'Tissue' column
    df['Tissue'] = [tissue1] * len(samples_tissue1) + [tissue2] * len(samples_tissue2)

    # Save the DataFrame to a new CSV file
    df.to_csv(f"../differential_expression_analysis/{tissue1}_{tissue2}_meta.txt", index=False)
    


def aggregation(gene_type, tissue1, tissue2, tissue1_df, tissue2_df, significant_results_file):
    """
    Aggregate gene regulatory network data for two tissues and differential expression analysis result.

    Args:
        gene_type (str): Type of gene ('TF_gene' or 'target_gene').
        tissue1 (str): Name of the first tissue.
        tissue2 (str): Name of the second tissue.
        tissue1_df (pd.DataFrame): DataFrame for the first tissue.
        tissue2_df (pd.DataFrame): DataFrame for the second tissue.
        significant_results_file (str): Path to the file with significant results.

    Returns:
        pd.DataFrame: Aggregated and merged data for both tissues.
    """
    # Sort and select top 3000 rows based on importance
    tissue1_df = tissue1_df.sort_values(by='importance', ascending=False).head(3000)
    tissue2_df = tissue2_df.sort_values(by='importance', ascending=False).head(3000)

    # Define the merge column based on the gene type
    merge_column = 'TF_gene' if gene_type == 'TF_gene' else 'target_gene'

    # Keep the first row after removing duplicates
    tissue1_df = tissue1_df.drop_duplicates(subset=merge_column, keep='first')
    tissue2_df = tissue2_df.drop_duplicates(subset=merge_column, keep='first')


    # Merge the tissue1 and tissue2 dataframes on the specified column
    common_df = pd.merge(tissue1_df, tissue2_df, on=merge_column, how='inner', suffixes=('_' + tissue1, '_' + tissue2))


    # Load differential expression analysis results
    diff_expression_df = pd.read_csv(significant_results_file)

    # Rename 'Unnamed: 0' column to match merge_column for merging
    diff_expression_df = diff_expression_df.rename(columns={'Unnamed: 0': 'TF_gene'})

    # Prepare the columns to select
    columns_to_select = ['TF_gene', 'log2FoldChange', 'regulation_status']

    if gene_type == 'TF_gene':
        # Determine which rows in diff_expression_df have TF_gene starting with 'ENST'
        enst_mask = diff_expression_df['TF_gene'].str.startswith('ENST')

        # Separate ENST and non-ENST rows
        enst_diff_df = diff_expression_df[enst_mask]
        non_enst_diff_df = diff_expression_df[~enst_mask]

        # Merge non-ENST rows based on TF_gene
        merged_df = pd.merge(common_df, non_enst_diff_df[columns_to_select], on='TF_gene', how='left')

        # Merge ENST rows based on TF_transcript_id_{tissue1} and TF_transcript_id_{tissue2}
        for tissue in [tissue1, tissue2]:
            enst_merged = pd.merge(common_df, enst_diff_df[columns_to_select],
                                   left_on=f'TF_transcript_id_{tissue}', right_on='TF_gene', how='left')
            # Handle cases with no match: drop 'TF_gene_x' and 'TF_gene_y' 
            if 'TF_gene_x' in enst_merged.columns:  # Check if 'TF_gene_x' exists (left column from merge)
                 enst_merged = enst_merged.drop('TF_gene_x', axis=1)  # Drop if no match
            if 'TF_gene_y' in enst_merged.columns:  # Check if 'TF_gene_y' exists (right column from merge)
                 enst_merged = enst_merged.drop('TF_gene_y', axis=1)  # Drop if no match
            else:
                 enst_merged = enst_merged.drop('TF_gene', axis=1)
            merged_df = merged_df.combine_first(enst_merged)
            # mark the flag as True to distinguish between TF based and Target based Merge
            merged_df['is_common_TF'] = True
    else:
        
        # Determine which rows in diff_expression_df have TF_gene starting with 'ENST'
        enst_mask = diff_expression_df['TF_gene'].str.startswith('ENST')

        # Separate ENST and non-ENST rows
        enst_diff_df = diff_expression_df[enst_mask]
        non_enst_diff_df = diff_expression_df[~enst_mask]

        # Merge non-ENST rows based on target_gene
        merged_df = pd.merge(common_df, non_enst_diff_df[columns_to_select], left_on='target_gene', right_on='TF_gene', how='left')
        merged_df = merged_df.drop('TF_gene', axis=1)

        # Merge ENST rows based on target_transcript_id_{tissue1} and target_transcript_id_{tissue2}
        for tissue in [tissue1, tissue2]:
            enst_merged = pd.merge(common_df, enst_diff_df[columns_to_select],
                                   left_on=f'target_transcript_id_{tissue}', right_on='TF_gene', how='left')
            enst_merged = enst_merged.drop('TF_gene', axis=1)
            merged_df = merged_df.combine_first(enst_merged)
            # mark the flag as False to distinguish between TF based and Target based Merge
            merged_df['is_common_TF'] = False
        
        # Rename columns if gene_type is 'target_gene'
        column_name_mapping = {
            'target_gene': 'TF_gene',
            f'target_transcript_id_{tissue1}': f'TF_transcript_id_{tissue1}',
            f'target_transcript_id_{tissue2}': f'TF_transcript_id_{tissue2}',
            f'TF_gene_{tissue1}': f'target_gene_{tissue1}',
            f'TF_transcript_id_{tissue1}': f'target_transcript_id_{tissue1}',
            f'TF_gene_{tissue2}': f'target_gene_{tissue2}',
            f'TF_transcript_id_{tissue2}': f'target_transcript_id_{tissue2}'
        }
        merged_df = merged_df.rename(columns=column_name_mapping)
        
    # Reorder columns
    new_columns_order = [
        'TF_gene', f'tissue_{tissue1}', f'TF_transcript_id_{tissue1}', f'target_gene_{tissue1}', f'target_transcript_id_{tissue1}', 
        f'importance_{tissue1}', f'tissue_{tissue2}', f'TF_transcript_id_{tissue2}', f'target_gene_{tissue2}', 
        f'target_transcript_id_{tissue2}', f'importance_{tissue2}', 'log2FoldChange', 'regulation_status','is_common_TF'
    ]
    merged_df = merged_df.reindex(columns=new_columns_order)
    merged_df.fillna(0, inplace=True)

    return merged_df


def map_names(df, gene_mapping_file, transcript_mapping_file, tissue1, tissue2):
    """
    Map gene and transcript IDs to their corresponding names.

    Args:
        df (pd.DataFrame): Input DataFrame with gene and transcript IDs.
        gene_mapping_file (str): Path to the gene mapping file.
        transcript_mapping_file (str): Path to the transcript mapping file.
        tissue1 (str): Name of the first tissue.
        tissue2 (str): Name of the second tissue.

    Returns:
        pd.DataFrame: DataFrame with mapped gene and transcript names.
    """
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
    gene_mapping_file = os.path.join(parent_dir, gene_mapping_file)
    transcript_mapping_file = os.path.join(parent_dir, transcript_mapping_file)
    
    # Load gene and transcript mapping datasets
    gene_mapping_df = pd.read_csv(gene_mapping_file)
    transcript_mapping_df = pd.read_csv(transcript_mapping_file)
    
    # Drop rows with NaN values
    gene_mapping_df.dropna(subset=['Gene name'], inplace=True)
    transcript_mapping_df.dropna(subset=['Transcript name'], inplace=True)
    
    # Create dictionaries for mapping
    gene_mapping_dict = dict(zip(gene_mapping_df['Gene stable ID'], gene_mapping_df['Gene name']))
    transcript_mapping_dict = dict(zip(transcript_mapping_df['Transcript stable ID'], transcript_mapping_df['Transcript name']))
    
    # Map gene stable IDs to gene names
    df['TF_gene'] = df['TF_gene'].map(gene_mapping_dict).fillna(df['TF_gene'])
    df[f'target_gene_{tissue1}'] = df[f'target_gene_{tissue1}'].map(gene_mapping_dict).fillna(df[f'target_gene_{tissue1}'])
    df[f'target_gene_{tissue2}'] = df[f'target_gene_{tissue2}'].map(gene_mapping_dict).fillna(df[f'target_gene_{tissue2}']) 
   
    # Function to map transcript IDs
    def map_transcript(transcript_id, tissue):
        if transcript_id != 0:
            return transcript_mapping_dict.get(transcript_id, transcript_id)
        else:
            return 0 if tissue not in transcript_mapping_dict.values() else transcript_id
    
    # Map transcript IDs to transcript names
    for tissue in [tissue1, tissue2]:
        df[f'TF_transcript_id_{tissue}'] = df[f'TF_transcript_id_{tissue}'].apply(lambda x: map_transcript(x, tissue))
        df[f'target_transcript_id_{tissue}'] = df[f'target_transcript_id_{tissue}'].apply(lambda x: map_transcript(x, tissue))
    
    # Output CSV file
    df.to_csv(f"../visualization/{tissue1}_{tissue2}_vis.csv", index=False)
    
    return df



def map_names_GSEA_target(gene_df, gene_mapping_file, transcript_mapping_file,diff_result):
    """
    Prepare input data for Gene Set Enrichment Analysis (GSEA).
    Select top 3000 target genes from GRNBoost2 network and combine with differential expression analysis results.
    Map Ensembl IDs to their external names.

    Args:
        gene_df (str): Path to the gene DataFrame file.
        gene_mapping_file (str): Path to the gene mapping file.
        transcript_mapping_file (str): Path to the transcript mapping file.
        diff_result (str): Path to the differential expression results file.

    Returns:
        pd.DataFrame: DataFrame prepared for GSEA with mapped gene names.
    """
    # Read GRN network
    df = pd.read_csv(gene_df)
    # select top 3000 based on importance
    df = df.sort_values(by='importance', ascending=False).head(3000)
    
    # Read differential expression analysis result
    de_results = pd.read_csv(diff_result)
    de_results = de_results.rename(columns={'Unnamed: 0': 'target_gene'})
    # Merge GRN and DGE based on target gene
    gene_list_df = pd.merge(df, de_results[['target_gene', 'log2FoldChange']], on="target_gene", how="inner")
    
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
    sys.path.append(parent_dir)
    gene_mapping_file = os.path.join(parent_dir, gene_mapping_file)
    transcript_mapping_file = os.path.join(parent_dir, transcript_mapping_file)
    
    # Load gene and transcript mapping datasets
    gene_mapping_df = pd.read_csv(gene_mapping_file)
    transcript_mapping_df = pd.read_csv(transcript_mapping_file)
    
    # Drop rows with NaN values
    gene_mapping_df.dropna(subset=['Gene name'], inplace=True)
    transcript_mapping_df.dropna(subset=['Transcript name'], inplace=True)
    
    # Create dictionaries for mapping
    gene_mapping_dict = dict(zip(gene_mapping_df['Gene stable ID'], gene_mapping_df['Gene name']))
    transcript_mapping_dict = dict(zip(transcript_mapping_df['Transcript stable ID'], transcript_mapping_df['Transcript name']))
    
    #If the ids start with ENSG, map gene-stable IDs to gene names
    #else map transcript-stable IDs to transcript names
    def map_gene_name(x):
        if x.startswith('ENSG'): 
            return gene_mapping_dict.get(x, 'Unknown')
        else:
            return transcript_mapping_dict.get(x, 'Unknown')

    gene_list_df['target_gene'] = gene_list_df['target_gene'].map(map_gene_name)
    gene_list_df.dropna(subset=['target_gene'], inplace=True)
    gene_list_df=gene_list_df.drop(['tissue', 'TF_gene','TF_transcript_id','target_transcript_id','importance'], axis=1)


    return gene_list_df


def run_gsea(tissue, diff_result, gene_set='GO_Biological_Process_2023', min_size=5, max_size=1000, permutation_num=1000, threads=4, seed=6):
    """
    Run Gene Set Enrichment Analysis (GSEA) for a specific tissue.
    This function performs GSEA to identify significantly enriched biological pathways in the gene regulatory network.

    Args:
        tissue (str): Name of the tissue.
        diff_result (str): Path to the differential expression results file.
        gene_set (str): Gene set to use for GSEA (default: 'GO_Biological_Process_2023').
        min_size (int): Minimum size of gene set to consider (default: 5).
        max_size (int): Maximum size of gene set to consider (default: 1000).
        permutation_num (int): Number of permutations for GSEA (default: 1000).
        threads (int): Number of threads to use (default: 4).
        seed (int): Random seed for reproducibility (default: 6).

    Returns:
        None. Saves GSEA results and plots.
    """
    # Read data
    tissue_df = f'{tissue}_df.csv'
    gene_mapping_file='../gene_mapping.csv'
    transcript_mapping_file='../transcript_mapping.csv'
    # Map transcript_mapping_file names using map_names_GSEA
    tissue_df_name = map_names_GSEA_target(tissue_df, gene_mapping_file, transcript_mapping_file, diff_result)

    # Save mapped names to files
    tissue_mapped_file = f"{tissue}_GSEA_target_gene.txt"
    tissue_df_name.to_csv(tissue_mapped_file, index=False, header=False)
    
    # Perform GSEA using gseapy
    tissue_rnk = pd.read_csv(tissue_mapped_file)
    pre_res = gp.prerank(rnk=tissue_rnk, 
                         gene_sets=gene_set,
                         threads=threads,
                         min_size=min_size,
                         max_size=max_size,
                         permutation_num=permutation_num,
                         outdir=None,  # Don't write to disk
                         seed=seed,
                         verbose=True)
    
    # Print top results
    print(pre_res.res2d.head(7))
    pre_res.res2d.to_csv(f"{tissue}_GSEA_result.txt")    
    # Plot the results
    terms = pre_res.res2d.Term
    
    axs = pre_res.plot(terms=terms[0:5],
                       show_ranking=True,
                       figsize=(3, 4))
    axs.figure.savefig(f"{tissue}_GSEA_dotplot.png", bbox_inches='tight')

    ax = dotplot(pre_res.res2d,
                 column="FDR q-val",
                 title=gene_set,
                 cmap=plt.cm.viridis,
                 size=6,
                 figsize=(4, 5), cutoff=0.05,  
                 show_ring=False)
    
    # Show the plot
    plt.show()



def isoform(tissue1, tissue2, expression_file_tissue1, expression_file_tissue2):
    """
    Analyze isoform expression patterns between two tissues.
    It calculates the percentage of each isoform's expression in each tissue, focusing on isoforms expressed in both tissues. 

    Args:
        tissue1 (str): Name of the first tissue.
        tissue2 (str): Name of the second tissue.
        expression_file_tissue1 (str): Path to the expression data file for tissue1.
        expression_file_tissue2 (str): Path to the expression data file for tissue2.

    Returns:
        pd.DataFrame: DataFrame containing isoform expression percentages for both tissues.
    """
    # Read expression data for tissue1
    tissue1_data = pd.read_csv(expression_file_tissue1, index_col=0, low_memory=False).T
    tissue1_data.reset_index(inplace=True)
    numerical_columns_tissue1 = tissue1_data.columns.difference(['index', 'Gene stable ID'])
    tissue1_data[numerical_columns_tissue1] = tissue1_data[numerical_columns_tissue1].apply(pd.to_numeric, errors='coerce')

    # Calculate total expression across numerical columns for tissue1
    total_tissue1 = tissue1_data[numerical_columns_tissue1].sum(axis=1, skipna=True)
    tissue1_data['Total'] = total_tissue1

    columns_to_drop_tissue1 = tissue1_data.columns.difference(['index', 'Gene stable ID', 'Total'])
    tissue1_data = tissue1_data.drop(columns_to_drop_tissue1, axis=1)
    tissue1_data['Tissue'] = tissue1

    # Read expression data for tissue2
    tissue2_data = pd.read_csv(expression_file_tissue2, index_col=0, low_memory=False).T
    tissue2_data.reset_index(inplace=True)
    numerical_columns_tissue2 = tissue2_data.columns.difference(['index', 'Gene stable ID'])
    tissue2_data[numerical_columns_tissue2] = tissue2_data[numerical_columns_tissue2].apply(pd.to_numeric, errors='coerce')

    # Calculate total expression across numerical columns for tissue2
    total_tissue2 = tissue2_data[numerical_columns_tissue2].sum(axis=1, skipna=True)
    tissue2_data['Total'] = total_tissue2

    columns_to_drop_tissue2 = tissue2_data.columns.difference(['index', 'Gene stable ID', 'Total'])
    tissue2_data = tissue2_data.drop(columns_to_drop_tissue2, axis=1)
    tissue2_data['Tissue'] = tissue2

    # Concatenate tissue1 and tissue2 dataframes
    combined_data = pd.concat([tissue1_data, tissue2_data], ignore_index=True)

    # Group by Gene stable ID and index, and calculate the number of unique tissues for each isoform
    grouped_data = combined_data.groupby(['Gene stable ID', 'index']).agg({
        'Tissue': 'nunique',
        'Total': 'sum'
    }).reset_index()

    # Filter out rows where an isoform is expressed in only one tissue
    multi_tissue_isoforms = grouped_data[grouped_data['Tissue'] > 1]

    # Merge back with the combined data to retain only rows with multi-tissue isoforms
    filtered_combined_data = pd.merge(combined_data, multi_tissue_isoforms[['Gene stable ID', 'index']], 
                                      on=['Gene stable ID', 'index'])

    # Calculate the percentage of each isoform's expression in each tissue
    filtered_combined_data['Percentage'] = filtered_combined_data.groupby(['Gene stable ID', 'index'])['Total'].transform(
        lambda x: x / x.sum()
    )

    # Prepare output dataframe
    output_data = filtered_combined_data[['Gene stable ID', 'index', 'Tissue', 'Percentage']]
    # Pivot the data to have separate columns for each tissue
  

    return output_data
    
def merge_grn_with_percentages(grn_table, grouped_table, tissue1, tissue2):
    """
    Merge GRN table with percentage information from grouped_table.
    
    Args:
    grn_table (pd.DataFrame): The original GRN table
    grouped_table (pd.DataFrame): The table with percentage information
    tissue1 (str): Name of the first tissue
    tissue2 (str): Name of the second tissue
    
    Returns:
    pd.DataFrame: Merged table with added percentage information
    """
    
    # Check if 'is_common_TF' column exists in grn_table
    if 'is_common_TF' not in grn_table.columns:
        raise ValueError("The 'is_common_TF' column is missing in the GRN table.")
    
    # Prepare grouped_table
    grouped_tissue1 = grouped_table[grouped_table['Tissue'] == tissue1].rename(columns={'Percentage': f'{tissue1}_percentage'})
    grouped_tissue2 = grouped_table[grouped_table['Tissue'] == tissue2].rename(columns={'Percentage': f'{tissue2}_percentage'})
    
    # Initialize an empty list to store the merged rows
    merged_rows = []
    
    # Iterate over each row in grn_table
    for _, row in grn_table.iterrows():
        is_common_TF = row['is_common_TF']
        
        # Merge based on TF_gene if is_common_TF is True, else merge based on target_gene
        if is_common_TF:
            merge_on_tissue1 = ['TF_gene', f'TF_transcript_id_{tissue1}']
            merge_on_tissue2 = ['TF_gene', f'TF_transcript_id_{tissue2}']
        else:
            merge_on_tissue1 = [f'target_gene_{tissue1}', f'target_transcript_id_{tissue1}']
            merge_on_tissue2 = [f'target_gene_{tissue2}', f'target_transcript_id_{tissue2}']
        
        # Merge with tissue1 data
        merged_row = pd.merge(row.to_frame().T, grouped_tissue1, 
                              left_on=merge_on_tissue1,
                              right_on=['Gene stable ID', 'index'],
                              how='left')
        
        # Merge with tissue2 data
        merged_row = pd.merge(merged_row, grouped_tissue2, 
                              left_on=merge_on_tissue2,
                              right_on=['Gene stable ID', 'index'],
                              how='left',
                              suffixes=('', '_tissue2'))
        
        # Calculate the other tissue percentages
        merged_row[f'{tissue1}_opposite_percentage'] = 1 - merged_row[f'{tissue1}_percentage']
        merged_row[f'{tissue2}_opposite_percentage'] = 1 - merged_row[f'{tissue2}_percentage']
        
        # Fill NaN values
        merged_row.fillna(0, inplace=True)
        
        # Append the merged row to the list
        merged_rows.append(merged_row)
    
    merged_table = pd.concat(merged_rows, ignore_index=True)
    
    # Reorder and select columns
    column_order = [
        'TF_gene', f'tissue_{tissue1}', f'TF_transcript_id_{tissue1}', f'target_gene_{tissue1}',
        f'target_transcript_id_{tissue1}', f'importance_{tissue1}', f'{tissue1}_percentage',
        f'{tissue1}_opposite_percentage', f'tissue_{tissue2}', f'TF_transcript_id_{tissue2}',
        f'target_gene_{tissue2}', f'target_transcript_id_{tissue2}', f'importance_{tissue2}',
        f'{tissue2}_percentage', f'{tissue2}_opposite_percentage',
        'log2FoldChange', 'regulation_status', 'is_common_TF'
    ]
    merged_table = merged_table[column_order]
    
    return merged_table


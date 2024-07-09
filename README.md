# Alternative splicing aware Gene Regulatory Network inference
This repository contains the code and data for analyzing alternative splicing-aware gene regulatory networks (GRNs) across different tissues. The analysis pipeline addresses the impact of alternative splicing on transcription factor (TF) isoforms, which is often overlooked in traditional GRN inference methods. By incorporating isoform expression data into GRN inference, this approach enhances our understanding of gene regulation and tissue-specific processes.
## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Data Download](#data-download)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Visualization](#visualization)
- [License](#license)
  
## An overview of the project
RNA-seq data from the GTEx project is preprocessed and analyzed using the GRNBoost2 algorithm, a tree-based machine learning algorithm selected for its accuracy and efficiency in GRN inference.
The inferred networks are combined with differential expression analysis to identify tissue-specific regulation by TF isoforms.
Gene Set Enrichment Analysis (GSEA) is performed to uncover overrepresented biological pathways and functions associated with the genes regulated by TF isoforms.
The approach is validated across tissue pairs, and the resulting GRNs are visualized using Cytoscape-cola.


## Requirements
- Python 3.9
- R (for differential expression analysis)
- Required Python libraries: pandas, arboreto, pybiomart, gseapy, matplotlib, numpy, dask, distributed
- Required R packages: DESeq2, ggplot2, reshape2
## installation
1. Clone this repository:
  ```
    git clone https://github.com/HyunjungWang/Thesis.git
    cd Thesis
  ```
2. Install required Python libraries:
  ```
   pip install pandas arboreto pybiomart gseapy matplotlib numpy dask distributed
  ```
3. Install required R packages:
  ```
   if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
   BiocManager::install("DESeq2")
   install.packages(c("ggplot2", "reshape2"))
  ```
## Data requirement
- Transcript expression data, gene expression data
  
  [!Transcript expression data structre](https://github.com/HyunjungWang/Thesis/assets/84681099/61f4b421-ab96-4022-b351-a1e6afc4c595)
  
  [!Gene expression data structure](https://github.com/HyunjungWang/Thesis/assets/84681099/2b5bf624-9216-4399-8540-91d29b3211ce)
- GTEX samples
- Meta data for differential expression analysis
  
  [!Meta data](https://github.com/HyunjungWang/Thesis/assets/84681099/0e5e2610-ff67-4e65-9b4b-a5ba5c13ef27)

## Data download
If you don't have required data structure, you can download data here.

https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

- Gene TPMs
  
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz

- Transcript TPMs
  
GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz

- GTEX samples
  
The sample data for all the tissue comparisons in this repository (brain-kidney, bladder-uterus, and cervix ectocervix-fallopian tube) are already included in their respective directories (brain_kidney/, bladder_uterus/, and cervix_fallopiantube/). You don't need to download any additional sample data for these tissue comparisons.
If you want to perform additional tissue comparisons in the future, you can download the corresponding sample data from the GTEx portal using the provided link.

Users can extract GTEX sample IDs with util.extract_and_save_gtex_columns()

## Repository Structure 
### util.py
Contains modular functions for data collection, preprocessing, GRN inference, alternative splicing-aware GRN construction, aggregation, GSEA, and isoform expression analysis.
### brain_kidney/
Folder containing the Jupyter Notebook and data files for the brain-kidney comparison analysis.
### bladder_uterus/
Folder containing the Jupyter Notebook and data files for the bladder-uterus comparison analysis.
### cervix_fallopiantube/ 
Folder containing the Jupyter Notebook and data files for the cervix ectocervix-fallopian tube comparison analysis.
### differential_expression_analysis/ 
Folder containing the R script and data files for differential gene expression analysis using DESeq2.

Note: Each tissue comparison folder contains a zipped file that needs to be unzipped before running the analysis.
The compressed file can be either in .zip or .7z format.
### visualization/
Folder containing the HTML file, cytoscape-cola-js, and CSV files for the GRN visualization tool.


## Usage
### Data Collection and Preprocessing:

Use read_transcript and read_gene functions to extract and organize raw data based on the selected columns provided in a separate data file..
Use preprocessing_data function to filter protein-coding genes and remove low expression genes.


### GRN Inference:

Use run_grnboost function to infer GRNs using the GRNBoost2 algorithm.


### Alternative Splicing-Aware GRN Construction:

Use grn_mapping and graph_data functions to create an alternative splicing-aware GRN structure.


### Differential Gene Expression Analysis:

Navigate to the differential_expression_analysis folder.
Run the R script to perform differential expression analysis between tissue pairs.

### Tissue-specific Isoform Expression Analysis:

Use isoform function to analyze isoform expression across tissues.

Use merge_grn_with_percentages function to combine isoform percentages with GRN.

### Gene Set Enrichment Analysis:

Use map_names_GSEA_target and run_gsea functions to perform GSEA.

### Visualization:
The GRN Visualization tool allows users to explore the alternative splicing-aware GRNs through a web browser. To use the visualization tool:
1. Ensure the GRN.html file and the visualization CSV files are in the visualization/ folder.
2. Open a command prompt, navigate to the visualization/ folder, start a local Python HTTP server:
  ```
    cd \Thesis\visualization
    python -m http.server
  ```
3. Open a web browser and go to http://localhost:8000/GRN.html
4. Select different tissue comparisons, adjust the input size, and choose between viewing only common TF genes or all data.





## Copyright Notice

This master's thesis, including all accompanying software, code, and data, is the intellectual property of [Hyunjung Wang], [2024]. All rights reserved. 

This work is submitted as part of the requirements for a master's degree at [Your University Name]. It is shared confidentially with the thesis advisor and evaluation committee. Any reproduction or use of the content of this thesis, in whole or in part, requires explicit permission from the author.

© [Hyunjung Wang] [2024]

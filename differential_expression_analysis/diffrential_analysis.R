# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# library(GenomeInfoDb)
library(DESeq2)
library(ggplot2)
library(reshape2)

# Function to perform analysis for a given option
perform_analysis <- function(input_file_path, sample_data_path, output_file_path) {
  # Read input files
  countData <- as.matrix(read.csv(input_file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE))

  colnames(countData) <- gsub("\\.", "-", colnames(countData))
  # Read the text file into a data frame
  sample_data<- as.matrix(read.csv(sample_data_path, header = TRUE))
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData=round(countData), colData=sample_data, design=~ Tissue)
  
  # Perform differential expression analysis
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  print(res)
  
  # Filter significant results
  log2FC_threshold <- 1.5
  pvalue_threshold <- 0.05
  
  complete_cases <- complete.cases(res$padj, res$log2FoldChange)
  significant_results <- res[complete_cases & res$padj < pvalue_threshold & abs(res$log2FoldChange) > log2FC_threshold, ]
  
  significant_results$regulation_status <- ifelse(significant_results$padj < pvalue_threshold, 
                                                  ifelse(significant_results$log2FoldChange > 0, "Upregulated", "Downregulated"), 
                                                  "Not Significant")
  
  # Print and save the significant results
  print(significant_results)
  write.csv(significant_results, file = output_file_path, row.names = TRUE)
  
  
  
  volcano_data <- as.data.frame(res)
  volcano_data$significant <- volcano_data$padj < pvalue_threshold & abs(volcano_data$log2FoldChange) > log2FC_threshold
  
  # Create the volcano plot
  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("black", "red"), labels = c("Not Significant", "Significant")) +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "blue") +
    labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-value)", color = "Significance") +
    ggtitle("Volcano Plot")
  
  # Print the volcano plot
  print(volcano_plot)
  ggsave(filename = output_volcano, plot = volcano_plot, width = 8, height = 6)  # Adjust width and height as needed
  

}

# Print available options
cat("Available options:\n")
cat("1. brain vs kidney\n")
cat("2. bladder vs uterus\n")
cat("3. cervix_ectocervix vs fallopian_tube\n")

# Prompt user to select an option
selected_option <- as.numeric(readline(prompt = "Enter the option number (1, 2, 3): "))

# Determine the selected option and call the corresponding analysis function
if (selected_option == 1) {
  input_file_path <- "brain_kidney_gene_expression_sum.csv"
  sample_data_path <- "brain_kidney_meta.txt"
  output_file_path <- "brain_kidney_significant_results.csv"
  output_volcano <- "brain_kidney_volcano.png"
  

} else if (selected_option == 2) {
  input_file_path <- "bladder_uterus_gene_expression_sum.csv"
  sample_data_path <- "bladder_uterus_meta.txt"
  output_file_path <- "bladder_uterus_significant_results.csv"
  output_volcano <- "bladder_uterus_volcano.png"
  

} else {
  input_file_path <- "cervix_ectocervix_fallopian_tube_gene_expression_sum.csv"
  sample_data_path <- "cervix_ectocervix_fallopian_tube_meta.txt"
  output_file_path <- "cervix_ectocervix_fallopian_tube_significant_results.csv"
  output_volcano <- "cervix_ectocervix_fallopian_tube_volcano.png"
  
}

# Perform analysis using the selected file paths
perform_analysis(input_file_path, sample_data_path, output_file_path)


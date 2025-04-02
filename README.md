# NSCLC Differential Gene Expression Analysis

## Overview

This repository contains the computational pipeline and results for analyzing differential gene expression (DEG) in **Non-Small Cell Lung Cancer (NSCLC)** using **RNA-seq data**. The main objective of this project is to identify and analyze the DEGs associated with NSCLC, focusing on how these genes vary based on age, smoking status, histology, tumor stage, and other clinical factors. The analysis was carried out using **DESeq2**, a widely-used method for RNA-seq differential expression analysis, and gene enrichment analysis tools.

This repository includes various R scripts, data files, and outputs from the analysis, as well as visualizations such as heatmaps for identifying gene expression patterns.

## Files and Directories

- **`Project.Functions.R`**: Custom R functions used in the analysis pipeline.
- **`Project.Main.R`**: Main R script used to run the differential gene expression analysis using DESeq2.
- **`Raw_Counts_GSE81089.tsv`**: Raw RNA-seq count data used for the DESeq2 analysis.
- **`FPKM_cufflinks.tsv`**: Alternative gene expression data in FPKM format.
- **`Significant_DEGs.tsv`**: Output file containing the results of the differential gene expression analysis, including the significant DEGs.
- **`Top_500_DEGs_by_log2foldchange.xlsx`**: Excel file listing the top 500 DEGs by log2 fold change.
- **`Sig_DEGS_Excel.xlsx`**: Merged file with significant DEGs data.
- **`MOET genes lung cancer source ensg.txt`**: List of genes identified as related to lung cancer from MOET analysis.
- **`SraRunTable Population data.xlsx`**: Patient and population data used in the analysis.
- **`moet_log2_EnsembleID.txt`**: Results from the MOET analysis for the top genes with log2 fold change.
- **`Heatmap analysis`**: R scripts and results for heatmap visualizations.
- **`.gitignore`**: Gitignore file to exclude unnecessary files from version control (e.g., .DS_Store).
- **`.Rhistory`**: R session history for tracking commands run during the analysis.

## Analysis Overview

1. **Data Preprocessing**: Raw RNA-seq count data was cleaned and processed to prepare it for differential expression analysis.
2. **Differential Expression Analysis**: DESeq2 was used to identify differentially expressed genes (DEGs) between NSCLC and healthy samples. Age, smoking status, histology, tumor stage, and other variables were considered to assess how these factors impact gene expression.
3. **Gene Enrichment Analysis**: Gene enrichment analysis was performed to identify associated pathways and diseases linked to the DEGs.
4. **Visualizations**: Heatmaps and other plots were generated to visualize gene expression patterns and clusters.
  
## Acknowledgments

We would like to thank Brian Bwanya, Maastricht University for providing the RNA-seq data and support for this project.


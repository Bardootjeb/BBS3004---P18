# Welcome to our project. For the code to work, a set of packages have to be
# installed. Follow the following commands to acquire the proper packages.

# Check if BiocManager is installed; install it if not
if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
} else {
  message("BiocManager is already installed")
}

# Create vector containing all the packages
packages <- c("DESeq2", "ggplot2", "dplyr", "pheatmap", 
"clusterProfiler", "org.Hs.eg.db", "GEOquery", "readr", "readxl")

# Check if required packages is installed; install it if not
if (!requireNamespace(packages)) {
  BiocManager::install(packages)
} else {
  message("Packages are already installed")
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Step 1. Preparatory Analysis

# Source files and data
source("Project.Functions.R")
require(DESeq2)
require(ggplot2)
require(dplyr)
require(pheatmap)
require(clusterProfiler)
require(org.Hs.eg.db)
require(GEOquery)
require(readr)  # For reading TSV files
require(readxl) # For reading Excel (XLSX) files

#load raw counts
counts <- read_tsv("GSE81089_raw_counts_GRCh38.p13_NCBI.tsv")
head(counts)

#load metadata
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

colnames(metadata1)

metadata.subset <- select(metadata, c(1,48, 49, 50, 51, 52, 53, 54, 56))
head(metadata.subset)

head(counts)
setdiff(colnames(counts), rownames(metadata.subset))  # Check if sample names in counts are in metadata

setdiff(rownames(metadata.subset), colnames(counts))  # Check if sample names in metadata are in counts

rownames(counts) <- counts[,1]  # Set first column (GeneID) as row names
counts <- counts[,-1]  # Remove the GeneID column from data

counts <- as.data.frame(counts)  # Convert tibble to a data frame
head(counts)  # Verify if the first column is GeneID
rownames(counts) <- counts[,1]  # Set the first column as row names
counts <- counts[,-1]  # Remove the GeneID column from the data
head(rownames(counts))  # Should display gene IDs
colnames(counts)  # Should display sample names

setdiff(colnames(counts), rownames(metadata))  # Should return character(0)
setdiff(rownames(metadata), colnames(counts))  # Should return character(0)



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
counts <- read.delim("FPKM_cufflinks.tsv", header=TRUE, 
                   row.names=1, sep="\t", check.names=FALSE)

head(counts)

# hoi

#load metadata
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
colnames(metadata) # I did colnames to see the different colomns

metadata.subset <- metadata[, c(1, 48, 49, 50, 51, 52, 53, 54, 56)]
head(metadata.subset)

rownames(metadata) <- metadata$title

colnames(counts)



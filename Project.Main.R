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
"clusterProfiler", "org.Hs.eg.db", "GEOquery")

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


#load raw counts
counts <- read.delim("FPKM_cufflinks.tsv", header=TRUE, 
                   row.names=1, sep="\t", check.names=FALSE)

# Looking at the head counts to see the type of data inside
head(counts)

#load metadata using the getGEO function
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)
#Extract metadata using pData function
metadata <- pData(phenoData(gse[[1]]))
# Look at the data inside. Head gives you the first 6
head(metadata)
# I did colnames to see the different colomns
colnames(metadata) 
# Create subset
metadata.subset <- metadata[, c(1, 48, 49, 50, 51, 52, 53, 54, 56)]
# What does this do??
rownames(metadata) <- metadata$title
# Look at the different names
colnames(counts)

# Extract the deaths
dead <- metadata.subset[, 3]
paste(dead)

# Ensure the output directory exists
if (!dir.exists("Output")) {
  dir.create("Output")
}


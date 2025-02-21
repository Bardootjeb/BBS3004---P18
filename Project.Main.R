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
metadata.subset <- metadata[, c(1, 8, 48, 49, 50, 51, 52, 53, 54, 56)]
# Look at the different names
colnames(counts)

# Renaming the colnames to the appropriate names to make it more readable
metadata.subset <- setNames(metadata.subset, c(
  "Title", "Source", "Age", "Life_Status", "Sex", "Histology", "Performance", 
  "Smoking_Status", "Tumor_stage", "Sample"
)[match(names(metadata.subset), c(
  "title", "source_name_ch1", "age:ch1", "dead:ch1", "gender:ch1", "histology:ch1", "ps who:ch1", 
  "smoking:ch1", "stage tnm:ch1", "tumor (t) or normal (n):ch1" 
  
))])

# Set column 'Sample' in metadata.subset as row names in metadata.subset (*to be able to match it later for deseq2 to column names of counts)
rownames(metadata.subset) <- metadata.subset$Sample

# Extract the deaths
dead <- metadata.subset[, 3]
paste(dead)

# Remove the last row from Data
head(counts)
dim(counts)
counts <- counts[-nrow(counts), ]

# Check if the last row is removed
dim(counts)  # Check new dimensions

# Ensure the output directory exists
if (!dir.exists("Output")) {
  dir.create("Output")
}

# Select our the genes of interest
interest.genes <- c("ENSG00000157764", "ENSG00000133703")

# Subset our genes of interest into new df by filtering on columns
express <- counts[rownames(counts) %in% interest.genes, , drop = FALSE]

# Convert the data into data frame
express <- as.data.frame(express)

# Reshape the expression data to better fit the dataframe
express$Gene <- rownames(express)
express<- reshape2::melt(express, id.vars = "Gene", variable.name = "Sample", 
                         value.name = "Expression")

# Merge with metadata
expression <- merge(express, metadata.subset, by = "Sample", all.x = TRUE)

# Plot expression levels of selected genes
save.pdf(function(){
  ggplot(expression, aes(x = Sample, y = Expression, fill = Gene)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Gene, scales = "free_y") +
    coord_flip() +  # Flip x and y axes
    theme_minimal() +
    labs(title = "Gene Expression Levels Across Samples",
         x = "Expression Level",
         y = "Sample")  # Swap x and y labels accordingly
}, "Sample Gene Expression Levels")


# Plot expression levels of selected genes
save.pdf(function(){
  ggplot(expression, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))# Rotate sample labels
}, "Sample Gene Expression Levels")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Step 2. Differential Gene Expression Analysis

# Making sure the row names in metadata.subset matches to column names in counts
all(colnames(counts) %in% rownames(metadata.subset))

# Find Columns in counts That Are Not in metadata.subset
setdiff(colnames(counts), rownames(metadata.subset))

# Remove everything after the underscore in counts
colnames(counts) <- sub("_.*", "", colnames(counts))
colnames(counts)    #check if it is removed

# Check again if row names in metadata.subset matches to column names in counts
all(colnames(counts) %in% rownames(metadata.subset))

# Check if they are in the same order
all(colnames(counts) == rownames(metadata.subset))

# Reorder metadata.subset rows to match the column order in counts
metadata.subset <- metadata.subset[match(colnames(counts), rownames(metadata.subset)), , drop = FALSE]

# Check if they now match
all(colnames(counts) == rownames(metadata.subset))

# Check the values in the counts
summary(counts)

# Convert all data values to Absolute values. (Non-negative)
Data <- abs(Data)

# Round values to integers
Data <- round(Data)

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = Data,
                              colData = metadata,
                              design = ~ Source)

dds

# Quality control
# Remove genes with low counts (choose one)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

keep2 <- rowMeans(counts(dds)) >=10
dds <- dds[keep2,]

dds


# Set the factor level




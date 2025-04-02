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

# Source files and require data
source("Project.Functions.R")
require(DESeq2)
require(ggplot2)
require(dplyr)
require(GEOquery)
require(tidyr)
require(pheatmap)

# Load FPKM normalized data 
FPKM_data <- read.delim("FPKM_cufflinks.tsv", header=TRUE, 
                        row.names=1, sep="\t", check.names=FALSE)

# Load metadata using the getGEO function
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)

# Extract metadata using pData function
metadata <- pData(phenoData(gse[[1]]))

# Create subset
metadata.subset <- metadata[, c(1, 8, 48, 49, 50, 51, 52, 53, 54, 56)]

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

# Remove the last row from FPKM_data
FPKM_data <- FPKM_data[-nrow(FPKM_data), ]

# Ensure the output directory exists
if (!dir.exists("Output")) {
  dir.create("Output")
}

# Select our the genes of interest B-RAF, K-RAS, EGFR
interest.genes <- c("ENSG00000157764", "ENSG00000133703", "ENSG00000146648")

# Subset our genes of interest into new df by filtering on columns
express <- FPKM_data[rownames(FPKM_data) %in% interest.genes, , drop = FALSE]

# Convert the data into data frame
express <- as.data.frame(express)

# Reshape the expression data to better fit the dataframe
express$Gene <- rownames(express)
express<- reshape2::melt(express, id.vars = "Gene", variable.name = "Sample", 
                         value.name = "Expression")

# Merge with metadata
express <- merge(express, metadata.subset, by = "Sample", all.x = TRUE)

# Plot expression levels of selected genes
# 3 genes separately in 3 plots
gene_colors <- c("pink", "lightblue", "lightgreen")  #create colour data
names(gene_colors) <- interest.genes  #assign a colour to each gene of interest

for(i in interest.genes){eplot <- ggplot(express %>% filter(Gene == i), aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge", fill = gene_colors[i]) +
  facet_wrap(~ Gene, scales = "free_y") +
  coord_flip() +  # Flip x and y axes
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples 2",
       x = "Expression Level",
       y = "Sample") 
print(eplot)}

# 3 genes together in one plot
save.pdf(function(){
  ggplot(express, aes(x = Sample, y = Expression, fill = Gene)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
    theme_minimal() +
    labs(title = "Gene Expression Levels Across Samples",
         x = "Sample",
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels
}, "Gene Expression")

# Boxplot of gene expression grouped by source 
save.pdf(function(){
  ggplot(express, aes(x = Source, y = Expression, fill = Source)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
    geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
    facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
    theme_minimal() +
    labs(title = "Gene Expression Levels by Source",
         x = "Source",
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability
}, "Boxplot; Source")

# Boxplot of gene expression grouped by sex
save.pdf(function(){
  ggplot(express, aes(x = Sex, y = Expression, fill = Sex)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
    geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
    facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
    theme_minimal() +
    labs(title = "Gene Expression Levels by Sex",
         x = "Sex",
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability
}, "Boxplot; Sex")

# Boxplot of gene expression grouped by smoking status
save.pdf(function(){
  ggplot(express, aes(x = Smoking_Status, y = Expression, fill = Smoking_Status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
    geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
    facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
    theme_minimal() +
    labs(title = "Gene Expression Levels by Smoking Status",
         x = "Smoking Status",
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability
}, "Boxplot; Smoking Status")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# Step 2. Differential Gene Expression Analysis

# Load the raw counts 
raw_counts <- read.delim("Raw_Counts_GSE81089.tsv", header=TRUE, 
                         row.names=1, sep="\t", check.names=FALSE)

# Making sure the row names in metadata.subset matches to column names in raw_counts
all(colnames(raw_counts) %in% rownames(metadata.subset))

# Reorder metadata.subset rows to match the column order in raw_counts
metadata.subset <- metadata.subset[match(colnames(raw_counts), rownames(metadata.subset)), , drop = FALSE]

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = metadata.subset,
                              design = ~ Source)

# Quality control
# Remove genes with low counts
keep <- rowMeans(counts(dds)) >=10
dds <- dds[keep,]

# Set the Source as factor instead of character
metadata.subset$Source <- as.factor(metadata.subset$Source)

# Sets the human non-malignant tissue as the base for when comparing
metadata.subset$Source <- relevel(metadata.subset$Source, ref = "Human non-malignant tissue")

# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)

# Extract results for Malignant Tissue vs. Human Tissue
res <- results(dds, contrast = c("Source", "Human NSCLC tissue", "Human non-malignant tissue" ))

# Filter for genes with padj < 0.01 (statistically significant) and log2FoldChange > 1 or < -1 (biologically meaningful)
deg_genes <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

# Check how many significant DEGs were found
nrow(deg_genes)

# Save results to a TSV file for further analysis
write.table(deg_genes, file= "Significant_DEGs.tsv", sep = "\t", col.names = F)

#making plots

# Convert results to a dataframe
res_df <- as.data.frame(res)

# Create a column for significance
res_df$significance <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1,
                              ifelse(res_df$log2FoldChange > 1, "Upregulated", "Downregulated"),
                              "Not Significant")

# Plot Volcano Plot
save.pdf(function(){
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +  # Vertical cutoff lines
    geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "black") +  # Horizontal cutoff line
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = "Volcano Plot of DEGs in NSCLC Tissue vs. Normal Tissue", 
         x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-Value") +
    theme(legend.title = element_blank())
}, "Volcano Plot")

# Convert results to a dataframe
res_df <- as.data.frame(res)

# Create a column for significance
res_df$significance <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1,
                              ifelse(res_df$log2FoldChange > 1, "Upregulated", "Downregulated"),
                              "Not Significant")

#convert into df
deg_genes_df <- as.data.frame(deg_genes)

# Sort by log2foldchange (column 2) in ascending order
sorted_deg_genes <- deg_genes_df[order(deg_genes_df[,2]), ]

# Take the first 500 genes
top_500_genes <- sorted_deg_genes[1:500, ]


# Add ENSEMBL IDs as the first column
top_500_genes_with_ids <- cbind(ENSEMBL_ID = rownames(top_500_genes), top_500_genes)

# Install and load the writexl package if needed
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

# Export to Excel
write_xlsx(top_500_genes_with_ids, path = "Top_500_DEGs_by_log2foldchange.xlsx")

# Plot MA plot
save.pdf(function(){
  plotMA(res)
}, "MA plot")


# Extract the top 10 upregulated and top 10 downregulated genes 
top_upregulated_NSCLC <- head(deg_genes[order(-deg_genes$log2FoldChange), ], 10)
top_downregulated_NSCLC <- head(deg_genes[order(-deg_genes$log2FoldChange), ], 10)

top_upregulated_NSCLC_df <- as.data.frame(top_upregulated_NSCLC)
top_downregulated_NSCLC_df <- as.data.frame(top_downregulated_NSCLC)

# upregulated and downregulated genes
upregulated_deg_genes <- subset(deg_genes, padj < 0.01 & log2FoldChange > 0)
downregulated_deg_genes <- subset(deg_genes, padj < 0.01 & log2FoldChange < 0)


# DESeq 2 for histology - Silke
# 1. change variables from characters to factors for 'histology' 

# Function to get the results for Histology
DSQ2("Histology", "Control", "metadata = metadata.subset")

metadata.subset$Histology <- as.factor(metadata.subset$Histology)

# 2. Construct a DESeqDataSet object
dds_histology <- DESeqDataSetFromMatrix(countData = raw_counts,
                                        colData = metadata.subset,
                                        design = ~ Histology)

# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_histology)) >=10
dds_histology <- dds_histology[keep,]

# 4. Set 'control' as the Reference
dds_histology$Histology <- relevel(dds_histology$Histology, ref = "Control")

# 5. Run DESeq2
dds_histology <- DESeq(dds_histology) 

# 6. Extract DEGs for groups: current smokers, ex smokers, never smokers 
res_control_vs_SC <- results(dds_histology, contrast = c ("Histology", "1", "Control"))
res_control_vs_AC <- results(dds_histology, contrast = c ("Histology", "2", "Control"))
res_control_vs_LC <- results(dds_histology, contrast = c ("Histology", "3", "Control"))

# 7. Extract DEGs within each comparison individually
DEGs_control_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
DEGs_control_vs_AC <- res_control_vs_AC[which(res_control_vs_AC$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
DEGs_control_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

DEGs_control_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res_control_vs_SC$log2FoldChange) > 1), ]
DEGs_control_vs_AC <- res_control_vs_AC[which(res_control_vs_AC$padj < 0.01 & abs(res_control_vs_AC$log2FoldChange) > 1), ]
DEGs_control_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res_control_vs_LC$log2FoldChange) > 1), ]

DEGs_control_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res_control_vs_SC$log2FoldChange) > 1), ]
DEGs_control_vs_AC <- res_control_vs_AC[which(res_control_vs_AC$padj < 0.01 & abs(res_control_vs_AC$log2FoldChange) > 1), ]
DEGs_control_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res_control_vs_LC$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_control_vs_SC, file= "DEGs_control_vs_SC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_AC, file= "DEGs_control_vs_AC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_LC, file= "DEGs_control_vs_LC.tsv", sep = "\t", col.names = F)

# 9. Making a plot 

# 8. Save to TSV
write.table(DEGs_control_vs_SC, file= "DEGs_control_vs_SC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_AC, file= "DEGs_control_vs_AC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_LC, file= "DEGs_control_vs_LC.tsv", sep = "\t", col.names = F)

# upregulated and downregulated genes
upregulated_control_vs_AC <- subset(DEGs_control_vs_AC, padj < 0.01 & log2FoldChange > 0)
downregulated_control_vs_AC <- subset(DEGs_control_vs_AC, padj < 0.01 & log2FoldChange < 0)

upregulated_control_vs_SC <- subset(DEGs_control_vs_SC, padj < 0.01 & log2FoldChange > 0)
downregulated_control_vs_SC <- subset(DEGs_control_vs_SC, padj < 0.01 & log2FoldChange < 0)

upregulated_control_vs_LC <- subset(DEGs_control_vs_LC, padj < 0.01 & log2FoldChange > 0)
downregulated_control_vs_LC <- subset(DEGs_control_vs_LC, padj < 0.01 & log2FoldChange < 0)

# top 5 upregulated and downregulated genes
top5_upregulated_degs_SC <- head(DEGs_control_vs_SC[order(-DEGs_control_vs_SC$log2FoldChange), ], 5)
top5_downregulated_degs_SC <- head(DEGs_control_vs_SC[order(-DEGs_control_vs_SC$log2FoldChange), ], 5)

top5_upregulated_degs_SC_df<- as.data.frame(top5_upregulated_degs_SC)
top5_downregulated_degs_SCdf <- as.data.frame(top5_downregulated_degs_SC)


top5_upregulated_degs_AC <- head(DEGs_control_vs_AC[order(-DEGs_control_vs_AC$log2FoldChange), ], 5)
top5_downregulated_degs_AC <- head(DEGs_control_vs_AC[order(-DEGs_control_vs_AC$log2FoldChange), ], 5)

top5_upregulated_degs_AC_df<- as.data.frame(top5_upregulated_degs_AC)
top5_downregulated_degs_ACdf <- as.data.frame(top5_downregulated_degs_AC)


top5_upregulated_degs_LC <- head(DEGs_control_vs_LC[order(-DEGs_control_vs_LC$log2FoldChange), ], 5)
top5_downregulated_degs_LC <- head(DEGs_control_vs_LC[order(-DEGs_control_vs_LC$log2FoldChange), ], 5)

top5_upregulated_degs_LC_df<- as.data.frame(top5_upregulated_degs_LC)
top5_downregulated_degs_LCdf <- as.data.frame(top5_downregulated_degs_LC)


# DEseq analysis with adenocarcinoma as reference

# Function to get the results for smoking
DSQ2("Histology", "Control")

metadata.subset$Histology <- as.factor(metadata.subset$Histology)

# 2. Construct a DESeqDataSet object
dds_histology <- DESeqDataSetFromMatrix(countData = raw_counts,
                                        colData = metadata.subset,
                                        design = ~ Histology)

# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_histology)) >=10
dds_histology <- dds_histology[keep,]

# 4. Set 'control' as the Reference
dds_histology$Histology <- relevel(dds_histology$Histology, ref = "2")

# 5. Run DESeq2
dds_histology <- DESeq(dds_histology) 

# 6. Extract DEGs for groups: current smokers, ex smokers, never smokers 
res_AC_vs_SC <- results(dds_histology, contrast = c ("Histology", "1", "2"))
res_AC_vs_LC <- results(dds_histology, contrast = c ("Histology", "3", "2"))

# 7. Extract DEGs within each comparison individually
DEGs_AC_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
DEGs_AC_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

DEGs_AC_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res_control_vs_SC$log2FoldChange) > 1), ]
DEGs_AC_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res_control_vs_LC$log2FoldChange) > 1), ]

DEGs_AC_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res_control_vs_SC$log2FoldChange) > 1), ]
DEGs_AC_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res_control_vs_LC$log2FoldChange) > 1), ]

# write a table for DEGs
write.table(DEGs_AC_vs_SC, file= "DEGs_AC_vs_SC.tsv", sep = "\t", col.names = F)
write.table(DEGs_AC_vs_LC, file= "DEGs_AC_vs_LC.tsv", sep = "\t", col.names = F)

# make a plot
plot_volcano(res_AC_vs_LC, "Volcano plot AC vs LC")
plot_volcano(res_AC_vs_SC, "Volcano plot AC vs SC")

# upregulated and downregulated genes
upregulated_AC_vs_LC <- subset(DEGs_AC_vs_LC, padj < 0.01 & log2FoldChange > 0)
downregulated_AC_vs_LC <- subset(DEGs_AC_vs_LC, padj < 0.01 & log2FoldChange < 0)

upregulated_AC_vs_SC <- subset(DEGs_AC_vs_SC, padj < 0.01 & log2FoldChange > 0)
downregulated_AC_vs_SC <- subset(DEGs_AC_vs_SC, padj < 0.01 & log2FoldChange < 0)

# top 5 up and downregulated genes
top5_upregulated_AC_vs_LC <- head(DEGs_AC_vs_LC[order(-DEGs_AC_vs_LC$log2FoldChange), ], 5)
top5_downregulated_AC_vsLC <- head(DEGs_AC_vs_LC[order(-DEGs_AC_vs_LC$log2FoldChange), ], 5)

top5_upregulated_ACLC_df<- as.data.frame(top5_upregulated_AC_vs_LC)
top5_downregulated_ACLC_df <- as.data.frame(top5_downregulated_AC_vsLC)


# Convert results to a dataframe
res_df <- as.data.frame(res)

# Create a column for significance
res_df$significance <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1,
                              ifelse(res_df$log2FoldChange > 1, "Upregulated", "Downregulated"),
                              "Not Significant")


# For functioncal enrichment analysis, create the data frames with the top 500 DEGs for every comparison
#Control vs AC
#convert into df
control_vs_AC_df <- as.data.frame(res_control_vs_AC)

# Sort by log2foldchange (column 2) in ascending order
sorteddegs_Control_vs_AC <- control_vs_AC_df[order(control_vs_AC_df[,2]), ]

# Take the first 500 genes
top500_Control_vs_AC <- sorteddegs_Control_vs_AC[1:500, ]

# Add ENSEMBL IDs as the first column
ids_control_vs_AC <- cbind(ENSEMBL_ID = rownames(top500_Control_vs_AC), top500_Control_vs_AC)

# Install and load the writexl package if needed
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)


# Export to Excel
write_xlsx(ids_control_vs_AC, path = "top500_control_vs_AC.xlsx")



#Control vs LC
control_vs_LC_df <- as.data.frame(res_control_vs_LC)

# Sort by log2foldchange (column 2) in ascending order
sorteddegs_Control_vs_LC <- control_vs_LC_df[order(control_vs_LC_df[,2]), ]

# Take the first 500 genes
top500_Control_vs_LC <- sorteddegs_Control_vs_LC[1:500, ]


# Add ENSEMBL IDs as the first column
ids_control_vs_LC <- cbind(ENSEMBL_ID = rownames(top500_Control_vs_LC), top500_Control_vs_LC)


# Export to Excel
write_xlsx(ids_control_vs_LC, path = "top500_control_vs_LC.xlsx")



#Control vs SC
control_vs_SC_df <- as.data.frame(res_control_vs_SC)

# Sort by log2foldchange (column 2) in ascending order
sorteddegs_Control_vs_SC <- control_vs_SC_df[order(control_vs_SC_df[,2]), ]

# Take the first 500 genes
top500_Control_vs_SC <- sorteddegs_Control_vs_SC[1:500, ]


# Add ENSEMBL IDs as the first column
ids_control_vs_SC <- cbind(ENSEMBL_ID = rownames(top500_Control_vs_SC), top500_Control_vs_SC)


# Export to Excel
write_xlsx(ids_control_vs_SC, path = "top500_control_vs_SC.xlsx")



#AC vs LC
AC_vs_LC_df <- as.data.frame(res_AC_vs_LC)

# Sort by log2foldchange (column 2) in ascending order
sorteddegs_AC_vs_LC <- AC_vs_LC_df[order(AC_vs_LC_df[,2]), ]

# Take the first 500 genes
top500_AC_vs_LC <- sorteddegs_AC_vs_LC[1:500, ]


# Add ENSEMBL IDs as the first column
ids_AC_vs_LC <- cbind(ENSEMBL_ID = rownames(top500_AC_vs_LC), top500_AC_vs_LC)


# Export to Excel
write_xlsx(ids_AC_vs_LC, path = "top500_AC_vs_LC.xlsx")


#AC vs SC
AC_vs_SC_df <- as.data.frame(res_AC_vs_SC)

# Sort by log2foldchange (column 2) in ascending order
sorteddegs_AC_vs_SC <- AC_vs_SC_df[order(AC_vs_SC_df[,2]), ]

# Take the first 500 genes
top500_AC_vs_SC <- sorteddegs_AC_vs_SC[1:500, ]


# Add ENSEMBL IDs as the first column
ids_AC_vs_SC <- cbind(ENSEMBL_ID = rownames(top500_AC_vs_SC), top500_AC_vs_SC)


# Export to Excel
write_xlsx(ids_AC_vs_SC, path = "top500_AC_vs_SC.xlsx")


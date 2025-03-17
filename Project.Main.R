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
require(tidyr)


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
expression <- merge(express, metadata.subset, by = "Sample", all.x = TRUE)

# Plot expression levels of selected genes
# 3 genes separately in 3 plots
gene_colors <- c("pink", "lightblue", "lightgreen")  #create colour data
names(gene_colors) <- interest.genes  #assign a colour to each gene of interest

for(i in interest.genes){eplot <- ggplot(expression %>% filter(Gene == i), aes(x = Sample, y = Expression, fill = Gene)) +
  geom_col(position = "dodge", fill = gene_colors[i]) +
  facet_wrap(~ Gene, scales = "free_y") +
  coord_flip() +  # Flip x and y axes
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples 2",
       x = "Expression Level",
       y = "Sample") 
print(eplot)}

# 3 genes together in one plot
ggplot(expression, aes(x = Sample, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Gene Expression Levels Across Samples",
       x = "Sample",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate sample labels

# Boxplot of gene expression grouped by source 
ggplot(expression, aes(x = Source, y = Expression, fill = Source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Source",
       x = "Source",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# Boxplot of gene expression grouped by sex
ggplot(expression, aes(x = Sex, y = Expression, fill = Sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Sex",
       x = "Sex",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# Boxplot of gene expression grouped by smoking status
ggplot(expression, aes(x = Smoking_Status, y = Expression, fill = Smoking_Status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot
  geom_jitter(width = 0.2, alpha = 0.6) +  # Adds individual points for visibility
  facet_wrap(~ Gene, scales = "free_y") +  # Separate plots for each gene
  theme_minimal() +
  labs(title = "Gene Expression Levels by Smoking Status",
       x = "Smoking Status",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# Step 2. Differential Gene Expression Analysis

# Load the raw counts 
raw_counts <- read.delim("Raw_Counts_GSE81089.tsv", header=TRUE, 
                     row.names=1, sep="\t", check.names=FALSE)

# Making sure the row names in metadata.subset matches to column names in raw_counts
all(colnames(raw_counts) %in% rownames(metadata.subset))

# Check if they are in the same order
all(colnames(raw_counts) == rownames(metadata.subset))

# Reorder metadata.subset rows to match the column order in raw_counts
metadata.subset <- metadata.subset[match(colnames(raw_counts), rownames(metadata.subset)), , drop = FALSE]

# Check if they now match
all(colnames(raw_counts) == rownames(metadata.subset))

# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = metadata.subset,
                              design = ~ Source)

print(dds)

# Quality control
# Remove genes with low counts
keep <- rowMeans(counts(dds)) >=10
dds <- dds[keep,]

print(dds)

# Set the Source as factor instead of character
metadata.subset$Source <- as.factor(metadata.subset$Source)

# Sets the human non-malignant tissue as the base for when comparing
metadata.subset$Source <- relevel(metadata.subset$Source, ref = "Human non-malignant tissue")

# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)

# Print a summary of DESeq2 results
print(dds)

# Extract results for Malignant Tissue vs. Human Tissue
res <- results(dds, contrast = c("Source", "Human non-malignant tissue", "Human NSCLC tissue" ))

# View a summary of the results
summary(res)

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
print(res_df$significance)

# Plot Volcano Plot
save.pdf(function(){
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
#  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme(legend.title = element_blank())
}, "Volcano Plot")


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DESeq2 for every variable

# First replace NA values with "Control" in metadata
metadata.subset<- metadata.subset%>%
  mutate(across(everything(), ~replace_na(.x, "Control"))) 

DSQ2(raw_counts, metadata.subset, "Smoking_Status", 3, "Smoking")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for smoking - Anne fleur
# 1. change variables from chracters to factors for 'smoking status' (or gender etc.) 
metadata.subset$Smoking_Status <- as.factor(metadata.subset$Smoking_Status)

# 2. Construct a DESeqDataSet object
dds_smoking <- DESeqDataSetFromMatrix(countData = raw_counts,
                                      colData = metadata.subset,
                                      design = ~ Smoking_Status)

# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_smoking)) >=10
dds_smoking <- dds_smoking[keep,]

# 4. Set never smokers (3) as the Reference
dds_smoking$Smoking_Status <- relevel(dds_smoking$Smoking_Status, ref = "3")

# 5. Run DESeq2
dds_smoking <- DESeq(dds_smoking)  

# 6. Extract DEGs for groups: current smokers, ex smokers, never smokers
res_never_vs_current <- results(dds_smoking, contrast = c ("Smoking_Status", "3", "1"))
res_never_vs_ex <- results(dds_smoking, contrast = c ("Smoking_Status", "3", "2"))
res_ex_vs_current <- results(dds_smoking, contrast = c ("Smoking_Status", "1", "2"))

# 7. Extract DEGs within each comparison individually
DEGs_never_vs_current <- res_never_vs_current[which(res_never_vs_current$padj < 0.01 & abs(res_never_vs_current$log2FoldChange) > 1), ]
DEGs_never_vs_ex <- res_never_vs_ex[which(res_never_vs_ex$padj < 0.01 & abs(res_never_vs_ex$log2FoldChange) > 1), ]
DEGs_ex_vs_current <- res_ex_vs_current[which(res_ex_vs_current$padj < 0.01 & abs(res_ex_vs_current$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_never_vs_current, file= "DEGs_never_vs_current.tsv", sep = "\t", col.names = F)
write.table(DEGs_never_vs_ex, file= "DEGs_never_vs_ex.tsv", sep = "\t", col.names = F)

# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq 2 for histology - Silke
# 1. change variables from characters to factors for 'histology' 
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
res_control_vs_SC <- results(dds_histology, contrast = c ("Histology", "Control", "1"))
res_control_vs_AC <- results(dds_histology, contrast = c ("Histology", "Control", "2"))
res_control_vs_LC <- results(dds_histology, contrast = c ("Histology", "Control", "3"))

# 7. Extract DEGs within each comparison individually
DEGs_control_vs_SC <- res_control_vs_SC[which(res_control_vs_SC$padj < 0.01 & abs(res_control_vs_SC$log2FoldChange) > 1), ]
DEGs_control_vs_AC <- res_control_vs_AC[which(res_control_vs_AC$padj < 0.01 & abs(res_control_vs_AC$log2FoldChange) > 1), ]
DEGs_control_vs_LC <- res_control_vs_LC[which(res_control_vs_LC$padj < 0.01 & abs(res_control_vs_LC$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_control_vs_SC, file= "DEGs_control_vs_SC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_AC, file= "DEGs_control_vs_AC.tsv", sep = "\t", col.names = F)
write.table(DEGs_control_vs_LC, file= "DEGs_control_vs_LC.tsv", sep = "\t", col.names = F)

# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for sex
# 1. change variables from chracters to factors for sex
metadata.subset$Sex <- as.factor(metadata.subset$Sex)

# 2. Construct a DESeqDataSet object
dds_sex <- DESeqDataSetFromMatrix(countData = raw_counts,
                                        colData = metadata.subset,
                                        design = ~ Sex)

# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_sex)) >=10
dds_sex <- dds_sex[keep,]

# 4. Set male as the Reference
dds_sex$Sex <- relevel(dds_sex$Sex, ref = "male")

# 5. Run DESeq2
dds_sex <- DESeq(dds_sex) 

# 6. Extract DEGs
res_male_vs_female <- results(dds_sex, contrast = c ("Sex", "male", "female"))

# 7. Extract DEGs within each comparison individually
DEGs_male_vs_female <- res_male_vs_female[which(res_male_vs_female$padj < 0.01 & abs(res_male_vs_female$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_male_vs_female, file= "DEGs_male_vs_female.tsv", sep = "\t", col.names = F)

# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for performance 
# 1. change variables from chracters to factors for performance
metadata.subset$Performance <- as.factor(metadata.subset$Performance)

# 2. Construct a DESeqDataSet object
dds_performance <- DESeqDataSetFromMatrix(countData = raw_counts,
                                  colData = metadata.subset,
                                  design = ~ Performance)
# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_performance)) >=10
dds_performance <- dds_performance[keep,]

# 4. Set '0 = normal activity' as the Reference
dds_performance$Performance <- relevel(dds_performance$Performance, ref = "0")

# 5. Run DESeq2
dds_performance <- DESeq(dds_performance)

# 6. Extract DEGs for groups:
# 7. Extract DEGs within each comparison individually
# 8. Save to TSV
# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for age
# 1. change variables from chracters to factors for performance 
metadata.subset$Age <- as.factor(metadata.subset$Age)

# 2. Construct a DESeqDataSet object
dds_age <- DESeqDataSetFromMatrix(countData = raw_counts,
                                         colData = metadata.subset,
                                         design = ~ Age)
# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_age)) >=10
dds_age <- dds_age[keep,]

# 4. Set '45 = youngest age' as the Reference ???

# 5. Run DESeq2
# 6. Extract DEGs for groups:
# 7. Extract DEGs within each comparison individually
# 8. Save to TSV
# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq 2 for tumor stage - Sabya
# 1. change variables from chracters to factors for tumor
metadata.subset$Tumor_stage <- as.factor(metadata.subset$Tumor_stage)

# 2. construct deseq set for Tumor_stage
dds_Tumorstage <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = metadata.subset,
                              design = ~ Tumor_stage)
# Quality control
# Remove genes with low counts
keep <- rowMeans(raw_counts(dds_Tumorstage)) >=10
dds_Tumorstage <- dds_Tumorstage[keep,]

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# Step 3. GO


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
res <- results(dds, contrast = c("Source", "Human non-malignant tissue", "Human NSCLC tissue" ))

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



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DESeq2 for every variable

# First replace NA values with "Control" in metadata
metadata.subset<- metadata.subset%>%
  mutate(across(everything(), ~replace_na(.x, "Control"))) 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for smoking - Anne fleur
# 1. change variables from chracters to factors for 'smoking status' (or gender etc.) 

# Function to get the results for smoking
DSQ2("Smoking_Status", 3, raw_counts, metadata.subset)

 ## The function makes the rest of the code obsolete
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
plot_volcano(res_ex_vs_current, "Volcano plot Ex vs Current Smokers")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq 2 for histology - Silke
# 1. change variables from characters to factors for 'histology' 

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
dds_histology$Histology <- relevel(dds_histology$Histology, ref = "Control")

# 5. Run DESeq2
dds_histology <- DESeq(dds_histology) 

# 6. Extract DEGs for groups: current smokers, ex smokers, never smokers 
res_control_vs_SC <- results(dds_histology, contrast = c ("Histology", "Control", "1"))
res_control_vs_AC <- results(dds_histology, contrast = c ("Histology", "Control", "2"))
res_control_vs_LC <- results(dds_histology, contrast = c ("Histology", "Control", "3"))

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
res_AC_vs_SC <- results(dds_histology, contrast = c ("Histology", "2", "1"))
res_AC_vs_LC <- results(dds_histology, contrast = c ("Histology", "2", "3"))

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
downregulated_AC_vs_SC <- subset(DEGs_AC_vs_LC, padj < 0.01 & log2FoldChange < 0)

upregulated_AC_vs_LC <- subset(DEGs_AC_vs_LC, padj < 0.01 & log2FoldChange > 0)
downregulated_AC_vs_SC <- subset(DEGs_AC_vs_LC, padj < 0.01 & log2FoldChange < 0)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for sex
# 1. change variables from chracters to factors for sex

# New function
DSQ2("Sex", "male")

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
DEGs_male_vs_female <- res_male_vs_female[which(res_male_vs_female$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_male_vs_female, file= "DEGs_male_vs_female.tsv", sep = "\t", col.names = F)

# 9. Making a plot 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for performance 

DSQ2("Performance", "0")

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

# Sort the dataset by Age
colData <- metadata.subset %>% arrange(Age)

# Assign variables as numeric
colData$Age <- as.numeric(colData$Age)

# Find the minimum and maximum ages
min_age <- min(colData$Age, na.rm = TRUE)  # Ensure NA values don't break it
max_age <- max(colData$Age, na.rm = TRUE)

# Filter for people within the youngest and oldest 10 year age range
youngest <- colData %>% filter(Age <= (min_age + 10))
oldest <- colData %>% filter(Age >= (max_age - 10))

# Combine the selected groups
selected_samples <- bind_rows(youngest, oldest)

# Convert Age into categorical variable
selected_samples$AgeGroup <- ifelse(selected_samples$Age <= (min_age + 10), "Young", "Old")
selected_samples$AgeGroup <- factor(selected_samples$AgeGroup, levels = c("Young", "Old"))

# Raw counts have to be adjusted
RCA <- raw_counts[, colnames(raw_counts) %in% selected_samples$Sample]

# Set row names of metadata to Sample IDs
rownames(selected_samples) <- selected_samples$Sample  

# Ensure RCA and metadata have the same order  
selected_samples <- selected_samples[colnames(RCA), ]

# function to analyse the rest automatically
DSQ2("AgeGroup", "Young", RCA, selected_samples)

# Subset count data to include only selected samples
RCA_subset <- raw_counts[, rownames(selected_samples)]

# Create DESeq2 dataset
dds_age <- DESeqDataSetFromMatrix(countData = RCA_subset,
                              colData = selected_samples,
                              design = ~ AgeGroup)

# 3. Quality control - Remove genes with low counts
keep <- rowMeans(counts(dds_age)) >=10
dds_age <- dds_age[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results comparing Old vs. Young
res <- results(dds, contrast = c("AgeGroup", "Old", "Young"))

# View top differentially expressed genes
resOrdered <- res[order(res$padj), ]
head(resOrdered)

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
keep <- rowMeans(counts(dds_Tumorstage)) >= 10
dds_Tumorstage <- dds_Tumorstage[keep, ]

# 4. Set 'control' as the Reference
dds_Tumorstage$Tumor_stage <- relevel(dds_Tumorstage$Tumor_stage, ref = "Control")


dds$Tumor_stage <- recode_factor(as.character(dds$Tumor_stage),
                                 "1" = "Stage 1",
                                 "2" = "Stage 1",
                                 "3" = "Stage 2",
                                 "4" = "Stage 2",
                                 "5" = "Stage 3",
                                 "6" = "Stage 3",
                                 "7" = "Stage 4")

dds$Tumor_stage <- factor(dds$Tumor_stage,
                          levels = c("Control", "Stage 1", "Stage 2", "Stage 3", "Stage 4"))

# 5. Run DESeq2
dds_Tumorstage <- DESeq(dds_Tumorstage) 

# 6. Extract DEGs for groups: current smokers, ex smokers, never smokers 
res_control_vs_stage1 <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "Control", "1"))
res_control_vs_stage2 <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "Control", "2"))
res_control_vs_stage3 <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "Control", "3"))
res_control_vs_stage4 <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "Control", "4"))

tumor_res <- c(res_control_vs_stage1, res_control_vs_stage2, res_control_vs_stage3, res_control_vs_stage4)

# 7. Extract DEGs within each comparison individually
DEGs_control_vs_stage1 <- res_control_vs_stage1[which(res_control_vs_stage1$padj < 0.01 & abs(res_control_vs_stage1$log2FoldChange) > 1), ]
DEGs_control_vs_stage2 <- res_control_vs_stage2[which(res_control_vs_stage2$padj < 0.01 & abs(res_control_vs_stage2$log2FoldChange) > 1), ]
DEGs_control_vs_stage3 <- res_control_vs_stage3[which(res_control_vs_stage3$padj < 0.01 & abs(res_control_vs_stage3$log2FoldChange) > 1), ]
DEGs_control_vs_stage4 <- res_control_vs_stage4[which(res_control_vs_stage4$padj < 0.01 & abs(res_control_vs_stage4$log2FoldChange) > 1), ]

# 8 Check how many significant DEGs were found
nrow(DEGs_control_vs_stage4)


# 9making plots

# Convert results to a dataframe
res_tumorstage <- as.data.frame(res_tumorstage)



# Create a column for significance
res_control_vs_stage1$significance <- ifelse(res_$padj < 0.01 & abs(res_df$log2FoldChange) > 1,
                                             ifelse(res_df$log2FoldChange > 1, "Upregulated", "Downregulated"),
                                             "Not Significant")
print(res_df$significance)

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

plot_volcano(res_control_vs_stage1, "tumorstage 1 vs control")
plot_volcano(res_control_vs_stage2, "tumorstage 2 vs control")
plot_volcano(res_control_vs_stage3, "tumorstage 3 vs control")
plot_volcano(res_control_vs_stage4, "tumorstage 4 vs control")
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
  # Step 3. GO

<<<<<<< HEAD
heatmap_plot(raw_counts, "sex", DEGs_male_vs_female)


# heatmap for 41 DEGs in NSCLC vs non-malignant, significantly associated with lung cancer 

# insert MOET output, ie the genes from the DEGs associated with lungdisease 
moetgenessource <- read.table("MOET genes lung cancer source ensg.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(moetgenessource) <- c("Gene_Symbol") 

# check overlap with expression data 
overlap_moetgenessource <- moetgenessource$Gene_Symbol[moetgenessource$Gene_Symbol %in% rownames(raw_counts)]

# Subset the data for MOET genes
expression_moetgenessource <- raw_counts[rownames(raw_counts) %in% overlap_moetgenessource, ]

# normalize it
log_overlap_moetgenessource <- log2(as.matrix(expression_moetgenessource) + 1)


# heatmap
pheatmap(log_overlap_moetgenessource, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row",
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "Heatmap of MOET Genes in Lung Cancer")



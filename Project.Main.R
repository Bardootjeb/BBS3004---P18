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
"clusterProfiler", "org.Hs.eg.db", "GEOquery", "openxlsx")

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
library(openxlsx)

# Load FPKM normalized data 
FPKM_data <- read.delim("FPKM_cufflinks.tsv", header=TRUE, 
                   row.names=1, sep="\t", check.names=FALSE)

# Load metadata using the getGEO function
gse <- getGEO(GEO = 'GSE81089', GSEMatrix = TRUE)

# Extract metadata using pData function
metadata <- pData(phenoData(gse[[1]]))

# Create subset
metadata.subset <- metadata[, c(1, 8, 48, 49, 50, 51, 52, 53, 54, 56)]

# Renaming the colnames of the metadata.subset to make it more readable
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

# Extract results for Human Tissue vs. Malignant Tissue
res <- results(dds, contrast = c("Source", "Human NSCLC tissue", "Human non-malignant tissue" ))

# Filter for genes with padj < 0.01 (statistically significant) and log2FoldChange > 1 or < -1 (biologically meaningful)
deg_genes <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

# Check how many significant DEGs were found
nrow(deg_genes)

# Save results to a TSV file for further analysis
write.table(deg_genes, file= "Significant_DEGs.tsv", sep = "\t", col.names = F)

#[making plot for the tumor vs normal tissue]

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DESeq2 for every variable

# First replace NA values with "Control" in metadata
metadata.subset<- metadata.subset%>%
  mutate(across(everything(), ~replace_na(.x, "Control"))) 

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for smoking - Anne fleur

# Function to get the results for smoking
DSQ2("Smoking_Status", 3, raw_counts, metadata.subset)

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
plot_volcano(res_ex_vs_current, "Volcano plot Ex vs Current Smokers")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq 2 for histology - Silke

# Function to get the results for histology
DSQ2("Histology", "Control", raw_counts, metadata.subset)

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


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for sex

# function for DESeq2 results for sex
DSQ2("Sex", "male", raw_counts, metadata.subset)

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
DEGs_male_vs_female <- res_male_vs_female[which(res_male_vs_female$padj < 0.01 & abs(res$log2FoldChange) > 1), ]

# 8. Save to TSV
write.table(DEGs_male_vs_female, file= "DEGs_male_vs_female.tsv", sep = "\t", col.names = F)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for performance 

#function to get results for performance
DSQ2("Performance", "0", raw_counts, metadata.subset)

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

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq2 for Age

# Sort the dataset by Age
colData <- metadata.subset %>% arrange(Age)

# Assign variables as numeric
colData$Age <- as.numeric(colData$Age)

# Find the minimum and maximum ages
min_age <- min(colData$Age, na.rm = TRUE)  # Ensure NA values don't break it
max_age <- max(colData$Age, na.rm = TRUE)

# Filter for people within the youngest and oldest 10 year age range
youngest <- colData %>% filter(Age <= (min_age + 13))
middle <- colData %>% filter(Age <= (min_age + 26))
oldest <- colData %>% filter(Age >= (max_age - 12))

# Combine the selected groups
selected_samples <- bind_rows(youngest, middle, oldest)

# Convert Age into categorical variable
selected_samples$AgeGroup <- case_when(
  selected_samples$Age <= (min(selected_samples$Age) + 13) ~ "Young",
  selected_samples$Age > (min(selected_samples$Age) + 13) & selected_samples$Age <= (max(selected_samples$Age) - 13) ~ "Middle",
  TRUE ~ "Old"  # Default case for remaining values
)

selected_samples <- selected_samples[!duplicated(selected_samples$Sample), ]

# Raw counts have to be adjusted
RCA <- raw_counts[, colnames(raw_counts) %in% selected_samples$Sample]

# Set row names of metadata to Sample IDs
rownames(selected_samples) <- selected_samples$Sample  

# Ensure RCA and metadata have the same order  
selected_samples <- selected_samples[colnames(RCA), ]

# function to analyse the rest automatically
DSQ2("AgeGroup", "Old", RCA, selected_samples)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# DESeq 2 for tumor stage - Sabya
# Step 1: Convert to character (just to avoid factor issues)
metadata.subset$Tumor_stage <- as.character(metadata.subset$Tumor_stage)

# Step 2: Recode sub-stages to single-digit stages
metadata.subset$Tumor_stage <- dplyr::recode(metadata.subset$Tumor_stage,
                                             "1" = "1", "2" = "1",
                                             "3" = "2", "4" = "2",
                                             "5" = "3", "6" = "3",
                                             "7"  = "4",
                                             "Control" = "Control"
)

# Step 3: Turn into factor, and set order if needed
metadata.subset$Tumor_stage <- factor(metadata.subset$Tumor_stage,
                                      levels = c("Control", "1", "2", "3", "4"))

# Step 4: Check it worked
table(metadata.subset$Tumor_stage)

# Step 5. construct deseq set for Tumor_stage
dds_Tumorstage <- DESeqDataSetFromMatrix(countData = raw_counts,
                                         colData = metadata.subset,
                                         design = ~ Tumor_stage)
# Step 6. Quality control
# Remove genes with low counts
keep <- rowMeans(counts(dds_Tumorstage)) >= 10
dds_Tumorstage <- dds_Tumorstage[keep, ]

# Step 7. Set 'control' as the Reference
dds_Tumorstage$Tumor_stage <- relevel(dds_Tumorstage$Tumor_stage, ref = "Control")


# Step 8. Run DESeq2
dds_Tumorstage <- DESeq(dds_Tumorstage) 

# Step 9. Extract DEGs for all stages
res_1_vs_control <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "1", "Control"))
res_2_vs_control <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "2", "Control"))
res_3_vs_control <- results(dds_Tumorstage, contrast = c ("Tumor_stage", "3", "Control"))


# Step 10. Extract DEGs within each comparison individually
DEGs_1_vs_control <- res_1_vs_control[which(res_1_vs_control$padj < 0.01 & abs(res_1_vs_control$log2FoldChange) > 1), ]
DEGs_2_vs_control <- res_2_vs_control[which(res_2_vs_control$padj < 0.01 & abs(res_2_vs_control$log2FoldChange) > 1), ]
DEGs_3_vs_control <- res_3_vs_control[which(res_3_vs_control$padj < 0.01 & abs(res_3_vs_control$log2FoldChange) > 1), ]


# Step 11. Check how many significant DEGs were found
nrow(DEGs_1_vs_control)
nrow(DEGs_2_vs_control)
nrow(DEGs_3_vs_control)


# step 12. Save DEG list
write.csv(DEGs_1_vs_control,
          file = "DEGs_1_vs_control.csv",
          row.names = TRUE)
write.csv(DEGs_2_vs_control,
          file = "DEGs_2_vs_control.csv",
          row.names = TRUE)
write.csv(DEGs_3_vs_control,
          file = "DEGs_3_vs_control.csv",
          row.names = TRUE)

# step 13. making plots with a function
plot_volcano(res_1_vs_control, "Stage 1 vs control")
plot_volcano(res_2_vs_control, "Stage 2 vs control")
plot_volcano(res_3_vs_control, "Stage 3 vs control")

#[comparison within tumorstages only]

# Subset metadata and counts to tumor stages only (1, 2, 3)
metadata_tumorOnly <- metadata.subset[metadata.subset$Tumor_stage %in% c("1", "2", "3"), ]
counts_tumorOnly <- raw_counts[, rownames(metadata_tumorOnly)]

# Make sure the stage is a factor
metadata_tumorOnly$Tumor_stage <- factor(metadata_tumorOnly$Tumor_stage,
                                         levels = c("1", "2", "3"))

# Create DESeqDataSet
dds_tumorOnly <- DESeqDataSetFromMatrix(
  countData = counts_tumorOnly,
  colData = metadata_tumorOnly,
  design = ~ Tumor_stage
)

# Filter low-expressed genes
dds_tumorOnly <- dds_tumorOnly[rowMeans(counts(dds_tumorOnly)) >= 10, ]

# Relevel to Stage 1
dds_stage1ref <- dds_tumorOnly
dds_stage1ref$Tumor_stage <- relevel(dds_stage1ref$Tumor_stage, ref = "1")
dds_stage1ref <- DESeq(dds_stage1ref)

# Extract results
res_2_vs_1 <- results(dds_stage1ref, contrast = c("Tumor_stage", "2", "1"))
res_3_vs_1 <- results(dds_stage1ref, contrast = c("Tumor_stage", "3", "1"))

# Filter DEGs
DEGs_2_vs_1 <- res_2_vs_1[which(res_2_vs_1$padj < 0.01 & abs(res_2_vs_1$log2FoldChange) > 1), ]
DEGs_3_vs_1 <- res_3_vs_1[which(res_3_vs_1$padj < 0.01 & abs(res_3_vs_1$log2FoldChange) > 1), ]

#Check how many significant DEGs were found
nrow(DEGs_2_vs_1)
nrow(DEGs_3_vs_1)


# Save DEG lists
write.table(rownames(DEGs_2_vs_1), "DEGs_Stage2_vs_Stage1.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(DEGs_3_vs_1), "DEGs_Stage3_vs_Stage1.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Relevel to Stage 2
dds_stage2ref <- dds_tumorOnly
dds_stage2ref$Tumor_stage <- relevel(dds_stage2ref$Tumor_stage, ref = "2")
dds_stage2ref <- DESeq(dds_stage2ref)

# Extract results
res_3_vs_2 <- results(dds_stage2ref, contrast = c("Tumor_stage", "3", "2"))

# Filter DEGs
DEGs_3_vs_2 <- res_3_vs_2[which(res_3_vs_2$padj < 0.01 & abs(res_3_vs_2$log2FoldChange) > 1), ]

# Check how many significant DEGs were found
nrow(DEGs_3_vs_2)

# Save DEG list
write.table(rownames(DEGs_3_vs_2), "DEGs_Stage3_vs_Stage2.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#making plots
plot_volcano(res_2_vs_1, "Stage 2 vs Stage 1")
plot_volcano(res_3_vs_1, "Stage 3 vs Stage 1")
plot_volcano(res_3_vs_2, "Stage 3 vs Stage 2")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
#HEALTHY TISSUES
#DESEQ2

# Becasue samples were only 19 and data was compromised (NA instead of actual data),we created a healhty data frame manually

healthy_metadata <- data.frame(
  sample = c("L511N", "L532N", "L561N", "L563N", "L566N", "L572N", "L606N", "L616N", "L644N", "L656N", "L661N",
             "L682N", "L723N", "L724N", "L736N",  "L738N",  "L809N", "L831N", "L881N"),  
  
  age = c(56, 69, 62,61, 58, 58, 55, 54, 67, 58, 66, 64, 67, 67, 67, 74, 72, 62, 70),
  
  Sex = c("female", "female", "female", "female", "male", "male", "male", "male", "male", 
          "male", "male", "female", "male", "female", "male", "female", "male", "female", "male")
)

rownames(healthy_metadata) <- healthy_metadata$sample

# Convert Age to numeric for classification
healthy_metadata$age <- as.numeric(as.character(healthy_metadata$age))

# Apply Age Group classification
healthy_metadata$Age <- case_when(
  healthy_metadata$age <= 58 ~ "Young",
  healthy_metadata$age >= 72 ~ "Old",
  TRUE ~ "Middle"
)

# Convert AgeGroup to a factor
healthy_metadata$age <- factor(healthy_metadata$Age, levels = c("Young", "Middle", "Old"))

# Filter raw_counts to include only the healthy samples
healthy_counts <- raw_counts[, rownames(healthy_metadata)]

# Check if sample matching worked
all(colnames(healthy_counts) == rownames(healthy_metadata))

DSQ2("Age", "Old", healthy_counts, healthy_metadata)



# Make sure the gender column is a factor
healthy_metadata$Sex <- as.factor(healthy_metadata$Sex)

# Create DESeq2 dataset for healthy samples only
dds_healthy_sex <- DESeqDataSetFromMatrix(
  countData = healthy_counts,
  colData = healthy_metadata,
  design = ~ Sex
)

# Filter low counts
keep <- rowMeans(counts(dds_healthy_sex)) >= 10
dds_healthy_sex <- dds_healthy_sex[keep, ]

# Set reference level
dds_healthy_sex$Sex <- relevel(dds_healthy_sex$Sex, ref = "female")  # Choose your reference

# Run DESeq
dds_healthy_sex <- DESeq(dds_healthy_sex)

# Extract results
res_healthy_sex <- results(dds_healthy_sex, contrast = c("Sex", "male", "female"))

# Extract significant DEGs
DEGs_healthy_sex <- res_healthy_sex[which(res_healthy_sex$padj < 0.01 & abs(res_healthy_sex$log2FoldChange) > 1), ]

# Convert to data frame and order by fold change
k_healthy <- as.data.frame(DEGs_healthy_sex)
orderedDEGs_healthy_sex <- k_healthy[order(k_healthy$log2FoldChange, decreasing = TRUE), ]

# View top DE genes
head(orderedDEGs_healthy_sex)

# Save results
write.table(DEGs_healthy_sex, file= "DEGs_healthy_sex.tsv", sep = "\t", col.names = F)
write_xlsx(DEGs_healthy_sex, path = "DEGs_healthy_sex.xlsx")

#comparison cancer and healthy genes 

# Get the gene lists from both datasets
healthy_sex_genes <- rownames(DEGs_healthy_sex)
cancer_sex_genes <- rownames(DEGs_male_vs_female)  
# Find overlapping genes
common_genes <- intersect(healthy_sex_genes, cancer_sex_genes)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
# Step 3. GO, heatmap of DEGs in NSCLC samples vs non malignant samples, sig. associated with ....

# Insert MOET output, ie the genes from the DEGs associated with lungdisease 
moetgenessource <- read.table("MOET genes lung cancer source ensg.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(moetgenessource) <- c("Gene_Symbol") 

# Check overlap with expression data 
overlap_moetgenessource <- moetgenessource$Gene_Symbol[moetgenessource$Gene_Symbol %in% rownames(raw_counts)]

# Subset the data for MOET genes
expression_moetgenessource <- raw_counts[rownames(raw_counts) %in% overlap_moetgenessource, ]

# Normalize it
log_overlap_moetgenessource <- log2(as.matrix(expression_moetgenessource) + 1)

# Split Heat Map into 2, Healthy & Non-Healthy
# Separate Tumor (“T”) and Normal (“N”) Samples
# Select columns ending with "T"
tumor_samples <- log_overlap_moetgenessource[, grep("T$", colnames(log_overlap_moetgenessource))]

# Select columns ending with "N"
normal_samples <- log_overlap_moetgenessource[, grep("N$", colnames(log_overlap_moetgenessource))]

### Put them All in One figure 
# Order columns: Normal (N) first, Tumor (T) second
ordered_columns <- c(grep("N$", colnames(log_overlap_moetgenessource), value = TRUE), 
                     grep("T$", colnames(log_overlap_moetgenessource), value = TRUE))
log_overlap_moetgenessource <- log_overlap_moetgenessource[, ordered_columns]

# Create a column annotation to distinguish Tumor vs. Normal
sample_types <- ifelse(grepl("T$", colnames(log_overlap_moetgenessource)), "Tumor", "Normal")
annotation_col <- data.frame(Type = factor(sample_types, levels = c("Normal", "Tumor")))
rownames(annotation_col) <- colnames(log_overlap_moetgenessource)

# Plot the heatmap
save.pdf(function(){
  pheatmap(log_overlap_moetgenessource, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           scale = "row",
           show_rownames = TRUE, 
           show_colnames = TRUE, 
           annotation_col = annotation_col,  # Highlight Tumor vs. Normal
           gaps_col = length(grep("N$", colnames(log_overlap_moetgenessource))), # Add gap between groups
           main = "Heatmap of MOET Genes (Normal vs. Tumor)")
}, "Heatmap of MOET Genes (Normal vs. Tumor)")


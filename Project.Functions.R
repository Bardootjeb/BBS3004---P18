# Function to save a plot to PDF. 
# plot_function must be a defined plotting function &
# filename will be the designated file name for the pdf
save.pdf <- function(plot_function, filename) {
  # Define the output PDF file path
  pdf_file <- file.path("Output", paste0(filename, ".pdf"))
  # Open the PDF device to save the plot
  pdf(pdf_file)
  # Call the plot function to create the plot
  print(plot_function())
  # Close the PDF device (this finalizes and writes the plot to the file)
  dev.off()
  # Inform the user that the plot has been saved
  message("Plot has been saved to: ", pdf_file)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Making a function to automate the deseq2 analysis
DSQ2 <- function(variable, ref_level, countdata, metadata){
  
  # Convert the variable to a factor
  metadata[[variable]] <- as.factor(metadata[[variable]])
  
  # Construct DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = countdata,
                                colData = metadata,
                                design = as.formula(paste("~", variable)))
  
  # Quality control - Remove genes with low counts
  keep <- rowMeans(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Set reference level
  dds[[variable]] <- relevel(dds[[variable]], ref = ref_level)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get levels for comparison
  levels_list <- levels(metadata[[variable]])
  
  # Create file to store results
  # Ensure the output directory exists
  if (!dir.exists(variable)) {
    dir.create(variable)
  }

  # Extract DEGs for all comparisons
  deg_results <- list()
  for (lvl in levels_list) {
    if (lvl != ref_level) {
      res <- results(dds, contrast = c(variable, ref_level, lvl))
      degs <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
      
      # Save results to file
      output_file <- file.path(variable, paste0(variable, "_", ref_level, "_vs_", lvl, ".tsv"))
      write.table(degs, file = output_file, sep = "\t", col.names = TRUE, row.names = TRUE)
      
      # Store results in a list
      deg_results[[paste0(ref_level, "_vs_", lvl)]] <- degs
      
      # Generate volcano plot
      plot_volcano(res, title = paste0(variable, "_", ref_level, "_vs_", lvl))
    }
  }
  
  return(deg_results)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Function to create a volcano plot of the degs

plot_volcano <- function(res, title) {
  save.pdf(function(){
  res_df <- as.data.frame(res)
  res_df$significance <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 1,
                                ifelse(res_df$log2FoldChange > 1, "Upregulated", "Downregulated"),
                                "Not Significant")
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black") +  # Vertical cutoff lines
    geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "black") +  # Horizontal cutoff line
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    theme(legend.title = element_blank())
  }, title)
}

# Function to create heatmaps
heatmap_plot <- function(counts_data, variable, DEGs){
  
# Subset expression data for DEGs related to your variable
heatmap_data <- counts_data[rownames(counts_data) %in% rownames(DEGs), ]

# Scale the data (row-wise normalization)
heatmap_data <- t(scale(t(heatmap_data)))

# Generate heatmap
pheatmap(heatmap_data, 
         annotation_col = metadata.subset[variable], 
         main = "Heatmap of DEGs", variable,
         color = colorRampPalette(c("blue", "white", "red"))(50))
}
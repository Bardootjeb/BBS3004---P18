# Function to save a plot to PDF. 
# plot_function must be a defined plotting function &
# filename will be the designated file name for the pdf
save.pdf <- function(plot_function, filename) {
  # Define the output PDF file path
  pdf_file <- file.path("Output", paste0(filename, ".pdf"))
  # Open the PDF device to save the plot
  pdf(pdf_file)
  # Call the plot function to create the plot
  plot_function()
  # Close the PDF device (this finalizes and writes the plot to the file)
  dev.off()
  # Inform the user that the plot has been saved
  message("Plot has been saved to: ", pdf_file)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
<<<<<<< HEAD
# Making a function to automate the deseq2 analysis

DSQ2 <- function(variable, ref_level) {
  
  # Convert the variable to a factor
  metadata.subset[[variable]] <- as.factor(metadata.subset[[variable]])
  
  # Construct DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                                colData = metadata.subset,
=======
# Making a function 

DSQ2 <- function(count_data, metadata, variable, ref_level, output_prefix) {
  
  # Convert the variable to a factor
  metadata[[variable]] <- as.factor(metadata[[variable]])
  
  # Construct DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = metadata,
>>>>>>> 97ef391a567d1fe6867099b56ffe048d27023cca
                                design = as.formula(paste("~", variable)))
  
  # Quality control - Remove genes with low counts
  keep <- rowMeans(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Set reference level
  dds[[variable]] <- relevel(dds[[variable]], ref = ref_level)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get levels for comparison
<<<<<<< HEAD
  levels_list <- levels(metadata.subset[[variable]])
  
  # Create file to store results
  # Ensure the output directory exists
  if (!dir.exists(variable)) {
    dir.create(variable)
  }
=======
  levels_list <- levels(metadata[[variable]])
>>>>>>> 97ef391a567d1fe6867099b56ffe048d27023cca
  
  # Extract DEGs for all comparisons
  deg_results <- list()
  for (lvl in levels_list) {
    if (lvl != ref_level) {
      res <- results(dds, contrast = c(variable, ref_level, lvl))
      degs <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1), ]
      
      # Save results to file
<<<<<<< HEAD
      output_file <- file.path(variable, paste0(variable, "_", ref_level, "_vs_", lvl, ".tsv"))
      write.table(degs, file = output_file, sep = "\t", col.names = TRUE, row.names = TRUE)
      
      # Store results in a list
      deg_results[[paste0(ref_level, "_vs_", lvl)]] <- degs
=======
      output_file <- paste0(output_prefix, "_", ref_level, "_vs_", lvl, ".tsv")
      write.table(degs, file = output_file, sep = "\t", col.names = TRUE, row.names = TRUE)
      
      # Store results in a list
      variable_results[[paste0(ref_level, "_vs_", lvl)]] <<- degs
>>>>>>> 97ef391a567d1fe6867099b56ffe048d27023cca
    }
  }
  
  return(deg_results)
}

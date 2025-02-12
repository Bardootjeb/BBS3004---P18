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

# Read the file
df <- read.delim("Significant_DEGs.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# Shuffle the rows of the dataframe
df <- df[sample(nrow(df)), ]

# Shuffle only the GeneExpression column (if needed)
if ("GeneExpression" %in% colnames(df)) {
  df$GeneExpression <- sample(df$GeneExpression)
}

# Save the modified dataframe
write.table(df, "Randomized_Significant_DEGs.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

# Print the first few rows to check
head(df)

# Put it in a file 
write.csv(df, "Randomized_Significant_DEGs.csv", row.names = TRUE)

head(df)

a <- 1:10
df$GeneExpression

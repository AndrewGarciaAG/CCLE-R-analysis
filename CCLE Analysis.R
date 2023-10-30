# Install the required packages if they are not already installed
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2")
}

# Set the working directory to a generic location
setwd("~/Documents")

# Read in the CSV files
gene_expression_data_1 <- read.csv("hncc.csv", header = TRUE)
gene_expression_data_2 <- read.csv("eso.csv", header = TRUE)
gene_expression_data_3 <- read.csv("lung.csv", header = TRUE)

# Remove the first 9 columns from each CSV file
gene_expression_data_1 <- gene_expression_data_1[, -c(1:9)]
gene_expression_data_2 <- gene_expression_data_2[, -c(1:9)]
gene_expression_data_3 <- gene_expression_data_3[, -c(1:9)]

# Convert each CSV file to a numeric matrix
gene_expression_data_1 <- as.matrix(gene_expression_data_1)
gene_expression_data_2 <- as.matrix(gene_expression_data_2)
gene_expression_data_3 <- as.matrix(gene_expression_data_3)

# Perform t-tests between each pair of gene expression matrices
he_t_test_results <- matrix(NA, ncol = 4, nrow = ncol(gene_expression_data_1))
colnames(he_t_test_results) <- c("t_statistic", "p_value", "df", "estimate")

for (i in 1:ncol(gene_expression_data_1)) {
  he_test_result <- t.test(gene_expression_data_1[, i], gene_expression_data_2[, i])
  he_t_test_results[i, "t_statistic"] <- he_test_result$statistic
  he_t_test_results[i, "p_value"] <- he_test_result$p.value
  he_t_test_results[i, "df"] <- he_test_result$parameter
  he_t_test_results[i, "estimate"] <- he_test_result$estimate[[1]] - he_test_result$estimate[[2]]
}

# Create a dot plot of the p-values from the t-tests
dot_plot <- ggplot(he_t_test_results, aes(x = colnames(he_t_test_results), y = p_value)) +
  geom_point(aes(color = colnames(he_t_test_results)), size = 5) +
  scale_y_log10() +
  theme_minimal() +
  labs(x = "Gene", y = "P-value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print the dot plot
print(dot_plot)
# Assuming you have a gene expression dataset in a variable called 'gene_set_data'

# Install and load the FactoMineR package
install.packages("FactoMineR")
library(FactoMineR)

gene_set_data <- read.csv("~/Documents/tempoSeq/gene counts csv/801_24h/801_24h_2nd_gene_counts.csv",
                          header = TRUE, row.names = 1)


gene_set_data <- data.frame(gene_set_data)
rownames(gene_set_data)
colnames(gene_set_data)
#delete all numbers behind the symbol "_"
gene_counts_isoform <- gsub( "_\\d+", " ", rownames(gene_set_data))
#generate a column named GeneIsoform
gene_set_data$GeneIsoform <- gene_counts_isoform


#count the duplicates of isoforms
sum(duplicated(gene_set_data$GeneIsoform))

colnames(gene_set_data)

#rename columns
colnames(gene_set_data)[1] <- "s1"
colnames(gene_set_data)[2] <- "s2"
colnames(gene_set_data)[3] <- "s3"
colnames(gene_set_data)[4] <- "s4"
colnames(gene_set_data)[5] <- "s5"
colnames(gene_set_data)[6] <- "s6"


#generate a form with combined isoforms 
gene_set_data1 <- aggregate(cbind(s1, s2, s3, s4, s5, s6) ~ GeneIsoform, gene_set_data, sum)

colnames(gene_set_data1)

# Assuming you have the 'gene_set_data1' data frame after the aggregation step

# Convert 'gene_set_data1' to a data frame with row names and column names
gene_set_data2 <- as.data.frame(gene_set_data1)

# Set row names based on the first column
rownames(gene_set_data2) <- gene_set_data1$GeneIsoform
gene_set_data2 <- gene_set_data2[,-1]

colnames(gene_set_data2)

# Identify rows with at least 3 values greater than 0 excluding the column of gene name
gene_set_data3 <- apply(gene_set_data2, 1, function(row) sum(row > 0) >= 3)
print(gene_set_data3)

# Subset the data frame to keep only the selected columns
gene_set_data4 <- gene_set_data2[gene_set_data3, ]


#save files
write_csv(gene_counts_cpm_merged_renamed_trimmed, 
          file = "Documents/tempoSeq/Tempo_Mar23_2023/Replicate_Apr20_2023/gene_counts_cpm_merged_renamed_trimmed.csv")

# Assuming your gene set data is stored in a matrix or dataframe called 'gene_set_data'

# Logarithmic transformation (base 2):when dealing with gene expression data measured on a logarithmic scale, such as microarray or RNA-seq data
transformed_data <- log2(gene_set_data4)

# Logarithmic transformation (base 10):when working with data on a decimal logarithmic scale
transformed_data <- log10(gene_set_data)

# Assuming your gene set data is stored in a variable called 'gene_set_data'

# Centering the data:Subtract the mean of each gene from the corresponding data points to center the data around zero.
centered_data <- scale(gene_set_data4, center = TRUE, scale = FALSE)

# Scaling the centered data:Divide the centered data by the standard deviation of each gene to scale the data and make genes comparable.
scaled_data <- scale(centered_data, center = FALSE, scale = TRUE)

# Assuming your gene set data is stored in a variable called 'gene_set_data'

# Remove rows with missing values
gene_set_data5 <- gene_set_data4[complete.cases(gene_set_data4), ]





library(caret)

# Assuming 'gene_set_data' contains the gene set data with a categorical variable 'category'

# Perform one-hot encoding
encoded_data <- predict(dummyVars(~ category, data = gene_set_data), newdata = gene_set_data)





# Perform PCA
pca_result <- PCA(gene_set_data, graph = FALSE)
pca_result <- PCA(gene_set_data5, graph = TRUE)

#1.Scree Plot
install.packages("factoextra")
library(factoextra)
fviz_eig(pca_result, addlabels = TRUE)
#2.Biplot:
fviz_pca_biplot(pca_result)
#3.Score Plot
plot(pca_result$ind$coord[, 1], pca_result$ind$coord[, 2], 
     xlab = "PC1", ylab = "PC2", main = "Score Plot")
#4.Loading Plot
plot(pca_result$var$coord[, 1], pca_result$var$coord[, 2], 
     xlab = "PC1", ylab = "PC2", main = "Loading Plot")
# Add labels to the plot
text(pca_result$var$coord[, 1], pca_result$var$coord[, 2], 
     labels = rownames(pca_result$var$coord), cex = 0.8)
#5.Scatter Plot of PC Scores
# Assuming you have performed PCA and obtained the PCA result in 'pca_result'

# Plot the scores for PC1 and PC2
plot(pca_result$ind$coord[, 1], pca_result$ind$coord[, 2], 
     xlab = "PC1", ylab = "PC2", main = "PCA Scatter Plot")

# Add labels to the data points (optional)
text(pca_result$ind$coord[, 1], pca_result$ind$coord[, 2], labels = rownames(pca_result$ind$coord), cex = 0.8)
#6.Biplot of PC Loadings and PC Scores
# Assuming you have performed PCA and obtained the PCA result in 'pca_result'

# Plot the biplot of PC loadings and PC scores
fviz_pca_biplot(pca_result, col.var = "blue", col.ind = "black", geom.ind = "point")
#7.Variance Explained Plot
# Assuming you have performed PCA and obtained the PCA result in 'pca_result'

# Plot the variance explained by each principal component
fviz_eig(pca_result, addlabels = TRUE)


# Access the results
# Eigenvalues
eigenvalues <- get_eigenvalue(pca_result)

# PC scores (coordinates of samples in the principal component space)
pc_scores <- get_pca_var(pca_result)$coord

# PC loadings (contribution of genes to each principal component)
pc_loadings <- get_pca_var(pca_result)$var
# Assuming you have performed PCA using FactoMineR and obtained the PCA result in 'pca_result'

# Access PC loadings
pc_loadings <- pca_result$ind$coord

# Continue with further analysis using the PC loadings

# Proportion of variance explained by each principal component
variance_explained <- get_eigenvalue(pca_result)$eig

# Assuming you have performed PCA using prcomp() and obtained the PCA result in 'pca_result'

# Access the eigenvalues for variance explained
variance_explained <- pca_result1$sdev^2 / sum(pca_result1$sdev^2)

# Print the variance explained by each principal component
print(variance_explained)


# Summary of the PCA results
summary(pca_result1)

# Assuming you have your gene expression data stored in a variable called 'gene_data'

# Perform PCA using prcomp()
pca_result1 <- prcomp(gene_set_data)

# Access PC loadings
pc_loadings <- pca_result1$x

# Continue with further analysis using the PC loadings


# Assuming you have performed PCA and obtained the PCA result in 'pca_result'

# Access PC loadings
pc_loadings <- get_pca_var(pca_result)$var
or
pc_loadings <- pca_result$rotation
or
pc_loadings <- pca_result$loadings
# Rank genes by loading magnitudes for a specific principal component (e.g., PC1)
pc1_loadings <- pc_loadings[, 1]  # Replace '1' with the desired PC index

# Sort genes by loading magnitudes for PC1 in descending order
sorted_genes <- names(sort(abs(pc1_loadings), decreasing = TRUE))

# Print the top genes contributing to PC1
num_top_genes <- 10  # Specify the number of top genes you want to analyze
top_genes <- sorted_genes[1:num_top_genes]
top_loadings <- pc1_loadings[top_genes]
for (i in 1:num_top_genes) {
  gene <- top_genes[i]
  loading <- top_loadings[i]
  print(paste("Gene:", gene, "| Loading:", loading))
}

# Visualize PC loadings
# Assuming you have a dataframe or matrix 'gene_expression' with gene expression data
# Plot loadings for the first two principal components
plot(pc_loadings[, 1], pc_loadings[, 2], xlab = "PC1 Loadings", ylab = "PC2 Loadings",
     xlim = c(-1, 1), ylim = c(-1, 1), main = "PC1 vs. PC2 Loadings")
text(pc_loadings[, 1], pc_loadings[, 2], labels = rownames(pc_loadings), cex = 0.8)

# Assuming you have performed PCA and obtained the PC loadings in 'pc_loadings'

# Specify the index of the principal component of interest (e.g., PC1)
pc_index <- 1

# Get the loadings for the specified principal component
pc_loadings <- pc_loadings[, pc_index]

# Calculate the absolute values of the loadings
abs_loadings <- abs(pc_loadings)

# Rank the genes by loading magnitudes in descending order
ranked_genes <- rownames(pc_loadings)[order(abs_loadings, decreasing = TRUE)]

# Print the ranked genes
print(ranked_genes)

dim(pc_loadings)

#Specify Directory
setwd("~/Data Landing Zone/Multiplex Data sets")
getwd()

#Load Libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(cluster)
library(caret)
library(tidyverse)
library(plot.matrix)
library(corrplot)
library(reshape2)
library(igraph)

#Load files
biomarkerFile <- "biomarkerThreshold.csv"
cellMeansFile <- "mean_cell_intensity.csv"
biomarkerThresholds <- read.csv (biomarkerFile)
cellMeans <- read.csv (cellMeansFile)

biomarkerThreshold <- data.frame(biomarkerThresholds)
cellMeans <- data.frame(cellMeans)

#extract biomarkers
unique_headers <- colnames(cellMeans)
biomarkers  <-
  unique_headers[-grep("parent", unique_headers, ignore.case = TRUE)]

#extract tissues
tissueList <- unique(cellMeans$Parent)

#remove rows with missing values
cellMeans <- na.omit(cellMeans)

# Filter data 
filtered <- cellMeans
for (biom in biomarkerThreshold$biomarker) {
  threshold <- biomarkerThreshold$threshold[biomarkerThreshold$biomarker == biom]
  filtered[, biom][cellMeans[, biom] < threshold] <- 0
  #filtered[, biom][cellMeans[, biom] >= threshold] <- 1
}

# PCA
pca_data <- prcomp(filtered[, biomarkers], scale = TRUE) 

# Clustering
d <- dist(pca_data$x)
hc <- hclust(d, method = "ward.D2")

# Plot dendrogram
plot(hc, hang = -1)

# Extract cluster assignments
clusters <- cutree(hc, k=5) 

# Add cluster to data
filtered$cluster <- clusters

pca_res <- prcomp(filtered[, biomarkers], scale = TRUE)

# Extract components
PC1 <- pca_res$x[,1] 
PC2 <- pca_res$x[,2]

# Add to data
filtered$PC1 <- PC1
filtered$PC2 <- PC2

# Plot
ggplot(filtered, aes(x=PC1, y=PC2, color=cluster)) +
  geom_point() +
  xlab("PC1") + 
  ylab("PC2") +
  ggtitle("PCA Cluster Plot")

#sub-set object such as to include 'all rows' and columns 2 to 8.
filteredValues <- filtered[ ,2:8] 



#Start here



corr_matrix <- cor(filteredValues)
binary_matrix <- as.matrix(filteredValues)
confusion_matrix <- confusionMatrix(binary_matrix)
plot(binary_matrix)

out <- do.call(rbind, lapply(list(binary_matrix), table))
rownames(out) <- c("Features")
barplot(t(out), beside = TRUE, legend.text = TRUE)


corrplot(cor(df))



# Replace upper triangle of the matrix with a placeholder value (e.g., 42)
my_cor_matrix[upper.tri(my_cor_matrix)] <- 42

# Melt the correlation matrix to create an adjacency list
my_cor_df <- melt(my_cor_matrix)
my_cor_df <- filter(my_cor_df, value != 42) %>% filter(Var1 != Var2)

# Rename columns for clarity
names(my_cor_df) <- c("from", "to", "weight")

# Print the first few rows of the adjacency list
head(my_cor_df)



# Create an igraph object from the adjacency list
net <- graph.data.frame(my_cor_df, directed = FALSE)

# Plot the graph (optional)
plot(net)
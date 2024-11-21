setwd("~/Data Landing Zone/SignalStar")
getwd()
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggdendro)


FILE234 <- "24in0465_individual_cluster_phenotype234.tsv"
data234 <- fread(FILE234) 

FILE91011 <- "24in0465_individual_cluster_phenotype91011.tsv"
data91011 <- fread(FILE91011) 

FILE131415 <- "24in0465_individual_cluster_phenotype131415.tsv"
data131415 <- fread(FILE131415) 

#join
combinedData <- rbind(data234, data91011, data131415)

#separate "name" by " " delimiter to make new columns
combinedData <- separate(combinedData, Name, c("study", "assay", "tissue", "plex", "section", "slide"), sep = " ")

#new dataframes, negative phenotype, coords, 4 relevent intensities
nucleiColumnNames <- c("section", grep("(Nuclei)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
nucleiData  <- na.omit(combinedData[, nucleiColumnNames])

#new dataframes, CD163 phenotype, coords, 4 relevent intensities
CD163ColumnNames <- c("section", grep("(CD163+)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
CD163Data <- na.omit(combinedData[, CD163ColumnNames])

#new dataframes, CD20 phenotype, coords, 4 relevent intensities
CD20ColumnNames <- c("section", grep("(CD20+)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
CD20Data <- na.omit(combinedData[, CD20ColumnNames])

#new dataframes, CD8 phenotype, coords, 4 relevent intensities
CD8ColumnNames  <- grep("CD8", names(combinedData), value = TRUE)

# exclude CD3 columns
cd8_not_cd3 <- CD8ColumnNames[!grepl("CD3", CD8ColumnNames)]

CD8ColumnNames <- c("section", cd8_not_cd3)

# Subset combinedData to new dataframe with selected columns
CD8Data <- na.omit(combinedData[, CD8ColumnNames])

#new dataframes, CD3 phenotype, coords, 4 relevent intensities
CD3ColumnNames  <- grep("CD3", names(combinedData), value = TRUE)

# exclude CD8 columns
cd3_not_cd8 <- CD3ColumnNames[!grepl("CD8", CD3ColumnNames)]

CD3ColumnNames <- c("section", cd3_not_cd8)

# Subset combinedData to new dataframe with selected columns
CD3Data <- na.omit(combinedData[, CD3ColumnNames])

#new dataframes, CD20 phenotype, coords, 4 relevent intensities
CD20ColumnNames <- c("section", grep("(CD20+)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
CD20Data <- na.omit(combinedData[, CD20ColumnNames])

#new dataframes, CD8+CD3+ phenotype, coords, 4 relevent intensities

CD8CD3ColumnNames <- c("section", grep("CD8\\+CD3\\+", names(combinedData), value=TRUE))

# Create new dataframe 
CD8CD3Data <- na.omit(combinedData[, CD8CD3ColumnNames])

#remove phenotypes from headers
colnames(CD8Data) <- gsub("\\s*\\(.*\\)", "", colnames(CD8Data))
colnames(CD8CD3Data) <- gsub("\\s*\\(.*\\)", "", colnames(CD8CD3Data))
colnames(CD3Data) <- gsub("\\s*\\(.*\\)", "", colnames(CD3Data))
colnames(CD20Data) <- gsub("\\s*\\(.*\\)", "", colnames(CD20Data))
colnames(CD163Data) <- gsub("\\s*\\(.*\\)", "", colnames(CD163Data))
colnames(nucleiData) <- gsub("\\s*\\(.*\\)", "", colnames(nucleiData))

#add phenotype column
CD8Data$phenotype <- 'CD8'
CD20Data$phenotype <- 'CD20'
nucleiData$phenotype <- 'Negative'
CD163Data$phenotype <- 'CD163'
CD3Data$phenotype <- 'CD3'
CD8CD3Data$phenotype <- 'CD8CD3'

#combine all individual cell data
cellData <- rbind(nucleiData, CD163Data, CD20Data, CD8Data, CD8CD3Data, CD3Data)

#extract list of sections
sectionList <- unique(cellData$section)

#Pick section to analyze----------------------------------------------------------------------------------------Input
nucleiSS <- nucleiData[nucleiData$section == "SS02", ]

allcellSS <- cellData[cellData$section == "SS02", ]

plot(nucleiSS$XCoord, nucleiSS$YCoord, 
     main="Original Data",
     xlab="XCoord", ylab="YCoord",
     pch=16)

# Perform k-means clustering
kmeans_models <- lapply(1:100, function(k){
  kmeans(nucleiSS[,c("XCoord", "YCoord")], k)
})

# Generate elbow plot 
wss <- sapply(kmeans_models, function(x) x$tot.withinss)
plot(1:100, wss,
     type="b", pch = 16, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Within groups sum of squares")

# Determine optimal number of clusters (from elbow plot)------------------
optimal_k <- 30

# Run k-means with optimal k 
final_model <- kmeans(nucleiSS[,c("XCoord", "YCoord")], optimal_k)

# Extract cluster assignments
clusters <- final_model$cluster

# Create color palette
cluster_colors <- rainbow(optimal_k)

# Plot data points colored by cluster 
plot(nucleiSS$XCoord, nucleiSS$YCoord, 
     col = cluster_colors[clusters],
     pch=16,
     main="Clustered Data",
     xlab="XCoord", ylab="YCoord")

# Add cluster centers 
points(final_model$centers[,1], final_model$centers[,2], 
       col=cluster_colors[clusters], pch=4, cex=2)

# Cluster centers 
centers <- final_model$centers

# Calculate distances 
dists <- matrix(NA, nrow=nrow(allcellSS), ncol=optimal_k)
for(i in 1:optimal_k){
  dists[,i] <- sqrt((allcellSS$XCoord - centers[i,1])^2 + 
                      (allcellSS$YCoord - centers[i,2])^2)
}

# Assign clusters
allcell_clusters <- apply(dists, 1, which.min)
allcellSS$cluster <- allcell_clusters

# Rainbow palette
cluster_colors <- rainbow(optimal_k) 

# Plot 
plot(allcellSS$XCoord, allcellSS$YCoord,
     col=cluster_colors[allcellSS$cluster], 
     pch=16,
     main="Clustered Data")

points(centers[,1], centers[,2],  
       col=cluster_colors, pch=4, cex=2)


#extract list of phenotypes
phenotypeList <- unique(allcellSS$phenotype)

#extract list of clusters
clusterList <- unique(allcellSS$cluster)

# Create an empty table to store the counts
totalCountTable <- matrix(0, nrow = length(clusterList), ncol = length(phenotypeList) + 1)

# Set the row names as cluster values
rownames(totalCountTable) <- clusterList

# Set the column names as phenotype values
colnames(totalCountTable) <- c(phenotypeList, "Total")

# Loop through each cluster and phenotype
for (i in 1:length(clusterList)) {
  for (j in 1:length(phenotypeList)) {
    # Count the number of rows that match the cluster and phenotype
    count <- sum(allcellSS$cluster == clusterList[i] & allcellSS$phenotype == phenotypeList[j])
    
    # Store the count in the table
    totalCountTable[i, j] <- count
  }
  
  # Calculate the total count for each cluster
  totalCountTable[i, length(phenotypeList) + 1] <- sum(totalCountTable[i, 1:length(phenotypeList)])
}

#make dataframe
totalCountTable <- as.data.frame(totalCountTable)

ryo <- colorRampPalette(c("yellow", "orange", "dark red"))(100)

# calculate ratios
ratioTable <- data.frame(apply(totalCountTable[,1:ncol(totalCountTable)-1], 2, 
                               function(x) x / totalCountTable$Total))

# remove negative
ratioTable <- ratioTable[ , colnames(ratioTable) != "Negative"]

# Hierarchical clustering
dendro <- hclust(dist(ratioTable), method = "complete")
hc_clusters <- cutree(dendro, k = 5)

# Melt dataframe 
melted <- melt(as.matrix(ratioTable))

# Plot
ggplot(melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = ryo) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Heatmap") +
  xlab("") +
  ylab("")

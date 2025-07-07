setwd("~/Data Landing Zone/SignalStar/25IN0124_clusters")
getwd()
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ggdendro)


#FILE <- "25IN0124_kidneys245.tsv"
FILE <- "25IN0124_livers2456.tsv"
#FILE <- "25IN0124_salivary24.tsv"

combinedData <- fread(FILE) 

#separate "name" by " " delimiter to make new columns
combinedData <- separate(combinedData, Name, c("study", "ID", "tissue", "RS", "Na", "originalstudy", "plex", "date", "slide"), sep = " ")

#new dataframes, negative phenotype, coords, ROI
nucleiColumnNames <- c("section", grep("(Nuclei)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
nucleiData  <- na.omit(combinedData[, nucleiColumnNames])

#new dataframes, CD163 phenotype, coords, ROI
CD163ColumnNames <- c("section", grep("(CD163+)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
CD163Data <- na.omit(combinedData[, CD163ColumnNames])

#new dataframes, CD20 phenotype, coords, ROI
CD20ColumnNames <- c("section", grep("(CD20+)", names(combinedData), value=TRUE))

# Subset combinedData to new dataframe with selected columns
CD20Data <- na.omit(combinedData[, CD20ColumnNames])

#new dataframes, CD8 phenotype, coords, ROI
CD8ColumnNames  <- grep("CD8", names(combinedData), value = TRUE)

# exclude CD3 columns
cd8_not_cd3 <- CD8ColumnNames[!grepl("CD3", CD8ColumnNames)]

CD8ColumnNames <- c("section", cd8_not_cd3)

# Subset combinedData to new dataframe with selected columns
CD8Data <- na.omit(combinedData[, CD8ColumnNames])

#new dataframes, CD3 phenotype, coords, ROI
CD3ColumnNames  <- grep("CD3", names(combinedData), value = TRUE)

# exclude CD8 columns
cd3_not_cd8 <- CD3ColumnNames[!grepl("CD8", CD3ColumnNames)]

CD3ColumnNames <- c("section", cd3_not_cd8)

# Subset combinedData to new dataframe with selected columns
CD3Data <- na.omit(combinedData[, CD3ColumnNames])


#new dataframes, CD8+CD3+ phenotype, coords, ROI

CD8CD3ColumnNames <- c("section", grep("CD8\\+CD3\\+", names(combinedData), value=TRUE))

# Create new dataframe 
CD8CD3Data <- na.omit(combinedData[, CD8CD3ColumnNames])

#remove phenotypes from headers
# Remove text before and including " - " from all column names
setnames(CD8Data, sub(".* - ", "", colnames(CD8Data)))
setnames(CD8CD3Data, sub(".* - ", "", colnames(CD8CD3Data)))
setnames(CD3Data,  sub(".* - ", "", colnames(CD3Data)))
setnames(CD20Data, sub(".* - ", "", colnames(CD20Data)))
setnames(CD163Data, sub(".* - ", "", colnames(CD163Data)))
setnames(nucleiData, sub(".* - ", "", colnames(nucleiData)))

#add phenotype column
CD8Data$phenotype <- 'CD8'
CD20Data$phenotype <- 'CD20'
nucleiData$phenotype <- 'Negative'
CD163Data$phenotype <- 'CD163'
CD3Data$phenotype <- 'CD3'
CD8CD3Data$phenotype <- 'CD8CD3'

#combine all individual cell data
cellData <- rbind(nucleiData, CD163Data, CD20Data, CD8Data, CD8CD3Data, CD3Data)


#extract list of section
sectionList <- unique(cellData$section)

#Pick section to analyze----------------------------------------------------------------------------------------Input
nucleiSS <- nucleiData[nucleiData$section == "SS04", ]

allcellSS <- cellData[cellData$section == "SS04", ]

plot(nucleiSS$"Envelope right", nucleiSS$"Envelope top",
     main="Original Data",
     xlab="X", ylab="Y")

# Get unique ROI names
ROIs <- unique(nucleiSS$ROIName)

# Create color palette
ROI_colors <- rainbow(length(ROIs))

# Map colors to the actual data points
color_map <- ROI_colors[as.factor(nucleiSS$ROIName)]

# Plot data points colored by ROI
plot(nucleiSS$"Envelope right", nucleiSS$"Envelope top",
     col = color_map,
     pch=16,
     main="ROIs",
     xlab="X", ylab="Y")


#extract list of phenotypes
phenotypeList <- unique(allcellSS$phenotype)

#extract list of ROIs
ROIList <- unique(allcellSS$ROIType)

# Create an empty table to store the counts
totalCountTable <- matrix(0, nrow = length(ROIList), ncol = length(phenotypeList) + 1)

# Set the row names as cluster values
rownames(totalCountTable) <- ROIList

# Set the column names as phenotype values
colnames(totalCountTable) <- c(phenotypeList, "Total")

# Loop through each ROI and phenotype
for (i in 1:length(ROIList)) {
  for (j in 1:length(phenotypeList)) {
    # Count the number of rows that match the ROI and phenotype
    count <- sum(allcellSS$ROIType == ROIList[i] & allcellSS$phenotype == phenotypeList[j])
    
    # Store the count in the table
    totalCountTable[i, j] <- count
  }
  
  # Calculate the total count for each cluster
  totalCountTable[i, length(phenotypeList) + 1] <- sum(totalCountTable[i, 1:length(phenotypeList)])
}

#make dataframe
totalCountTable <- as.data.frame(totalCountTable)

# Remove sum=1 rows
totalCountTable <- totalCountTable[rowSums(totalCountTable) != 2, ]

# calculate ratios
ratioTable <- data.frame(apply(totalCountTable[,1:ncol(totalCountTable)-1], 2, 
                               function(x) x / totalCountTable$Total))


# remove negative
ratioTable <- ratioTable[ , colnames(ratioTable) != "Negative"]

# Compute the distance matrix
dist_matrix <- dist(ratioTable, method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

# Cut the dendrogram into a desired number of clusters (e.g., 3)
clusters <- cutree(hc, k = 3)


# Plot the dendrogram
plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "", cex = 0.9)

# Melt dataframe 
melted <- melt(as.matrix(ratioTable))

# Convert 'value' column to numeric
melted$value <- as.numeric(melted$value)

# Specify colors
ryo <- colorRampPalette(c("yellow", "orange", "dark red"))(100)

# Plot Heatmap
ggplot(melted, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = ryo) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle("Heatmap") +
  xlab("") +
  ylab("")

# Add cluster assignments to the dataframe
ratioTable$cluster <- as.factor(clusters)


#New datatable
SSclusters <-allcellSS

# Extract row names and the "cluster" column from ratioTable
ratioTable_extracted <- data.frame(ROIType = rownames(ratioTable), cluster = ratioTable$cluster)

# Merge the extracted data with SSclusters based on the "ROIType" column
SSclusters <- merge(SSclusters, ratioTable_extracted, by = "ROIType", all.x = TRUE)


# Plot data points colored by cluster assignment
plot(SSclusters$"Envelope right", SSclusters$"Envelope top",
     col = SSclusters$cluster,
     pch=16,
     main="Hierarchical Clustering",
     xlab="X", ylab="Y")
                
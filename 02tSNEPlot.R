#Specify Directory
setwd("~/Data Landing Zone/Multiplex Data sets")
getwd()

#Load Libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Rtsne)

#Load files
biomarkerFile <- "biomarkerThreshold.csv"
cellMeansFile <- "mean_cell_intensity.csv"
biomarkerThresholds <- read.csv (biomarkerFile,  sep=',', header = T)
cellMeans <- read.csv (cellMeansFile, sep=',', header = T)

#format as dataframe
biomarkerThreshold <- data.frame(biomarkerThresholds)
cellMeans <- data.frame(cellMeans)

#remove rows with missing values
cellMeans <- na.omit(cellMeans)

#sub-set object such as to include 'all rows' and columns 2 to 8.
biomarkerValues <- cellMeans[ ,2:8] 
#sub-set object such as to include 'all rows' and column 1 parent 
biomarkerParent <- cellMeans[ ,1] 

# Create a color palette that matches the number of unique categories in biomarkerParent
unique_categories <- unique(biomarkerParent)
color_palette <- rainbow(length(unique_categories))

# Create a named vector of colors to match categories with colors
colors_for_plot <- setNames(color_palette, unique_categories)

tsne_results <- Rtsne(biomarkerValues, perplexity=30, check_duplicates = FALSE)

# Plotting
par(mfrow=c(1,2)) # To plot two images side-by-side
plot(tsne_results$Y, col = "blue", pch = 19, cex = 0.3)

# Second plot with points colored by biomarkerParent categories
plot(tsne_results$Y, col = colors_for_plot[biomarkerParent], pch = 19, cex = 0.3)

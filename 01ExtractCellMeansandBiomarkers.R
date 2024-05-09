setwd("~/Data Landing Zone/Multiplex Data sets")
getwd()
library(data.table)
library(dplyr)

#FILENAME <- "USL-2023-46567_1_Scene1_R1.ome.tif - USL-2023-46567_1_Scene1_R1.afi.txt"
#FILENAME <- "USL-2023-46567_1_Scene3_R1.ome.tif - USL-2023-46567_1_Scene3_R1.afi.txt"
FILENAME <- "USL-2023-46565_1.ome.tif - USL-2023-46565_1.afi.txt"

data <- read.table (FILENAME,  sep = "\t", header = TRUE)

orginal_headers <- names(data)

ignore_headers <-c("RED", "GREEN", "BLUE", "DAPI")

channels_to_ignore <- sapply(list(ignore_headers), paste, collapse="|")

#select relevant means and parent columns
selected_columns <-
  data %>% select(matches("mean|parent|object.ID", ignore.case = TRUE))
selected_columns <-
  selected_columns %>% select(-matches(channels_to_ignore))

cell_means <-
  selected_columns %>% select(matches("cell|parent", ignore.case = TRUE))

colnames(cell_means) <- 
  sub("\\.\\..*", "", colnames(cell_means))

#extract list of unique column names(biomarkers), excluding parent
unique_headers <- colnames(cell_means)
biomarkers  <-
  unique_headers[-grep("parent", unique_headers, ignore.case = TRUE)]

# Create a data frame
biomarkerThreshold <- data.frame(biomarker = biomarkers, threshold = rep(0, length(biomarkers)))

# Write the data frame to a CSV file
write.csv(biomarkerThreshold, file = "biomarkerThreshold.csv", row.names = FALSE)

# Define filename
additional_words <- "mean_cell_intensity"

# Create the filename
filename <- paste0( additional_words, ".csv")

# Write the dataframe to the CSV file
write.csv(cell_means, filename, row.names = FALSE)

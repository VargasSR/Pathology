#Specify Directory
setwd("~/Data Landing Zone/Raw_Spreadsheets")
getwd()

#Load Libraries
library(readxl)
library(data.table)
library(dplyr)
library(tidyverse)
library(writexl)

#list files
files <- list.files(path="~/Data Landing Zone/Raw_Spreadsheets", pattern="*.xls", full.names=TRUE)

# Initialize empty dataframe for SLIDE
slide_tab <- data.frame() 

# Loop through files
for(file in files) {
  
  # Get study name from filename
  study_name <- tools::file_path_sans_ext(basename(file))  
  
  # Read Slide data
  slide_data <- read_excel(file, sheet="SLIDE")
  
  # Modify SUBJID column
  slide_data$SUBJID <- paste0(slide_data$SUBJID, "-", study_name)
  
  # Add Study column
  slide_data$Study <- study_name
  
  # Bind to dataframe
  slide_tab <- rbind(slide_tab, slide_data)
  
}

#specify valid blocks and remove others
valid_blocks <- c("A1", "B1", "D1", "E1", "F1", "F2", "G2")
block_filter_slide_tab <- slide_tab[slide_tab$BLOCKID %in% valid_blocks,]

# Initialize empty dataframe for SUBJECT
subject_tab <- data.frame() 

# Loop through files
for(file in files) {
  
  # Get study name from filename
  study_name <- tools::file_path_sans_ext(basename(file))  
  
  # Read Slide data
  subject_data <- read_excel(file, sheet="SUBJECT")
  
  
  # Modify SUBJID column
  subject_data$SUBJID <- paste0(subject_data$SUBJID, "-", study_name)
  
  # Bind to dataframe
  subject_tab <- rbind(subject_tab, subject_data)
  
}

# Initialize empty dataframe for ORGAN
organ_tab <- data.frame() 

# Loop through files
for(file in files) {
  
  # Get study name from filename
  study_name <- tools::file_path_sans_ext(basename(file))  
  
  # Read Slide data
  organ_data <- read_excel(file, sheet="ORGAN")
  
  
  # Modify SUBJID column
  organ_data$SUBJID <- paste0(organ_data$SUBJID, "-", study_name)
  
  # Bind to dataframe
  organ_tab <- rbind(organ_tab, organ_data)
  
}

#specify valid organs and remove others
valid_organs <- c("Kidney", "Liver", "Lung", "Heart", "Gland, Thyroid", "Gland, Adrenal", "Thymus")
organ_filter_organ_tab <- organ_tab[organ_tab$MISPEC %in% valid_organs,]

# Initialize empty dataframe for MI FNDING
finding_tab <- data.frame() 

# Loop through files
for(file in files) {
  
  # Get study name from filename
  study_name <- tools::file_path_sans_ext(basename(file))  
  
  # Read Slide data
  finding_data <- read_excel(file, sheet="MI FINDING")
  
  # Add Study column
  finding_data$Study <- study_name
  
  # Modify SUBJID column
  finding_data$SUBJID <- paste0(finding_data$SUBJID, "-", study_name)
  
  # Bind to dataframe
  finding_tab <- rbind(finding_tab, finding_data)
  
}

#remove column
finding_tab <- finding_tab[, !names(finding_tab) %in% "Treatment Related Flag"]

# remove others invalid organs
organ_filter_finding_tab <- finding_tab[finding_tab$MISPEC %in% valid_organs,]


# Create a list of dataframes
tab_list <- list(block_filter_slide_tab, subject_tab, organ_filter_organ_tab, organ_filter_finding_tab) 

# Name each dataframe tab
names(tab_list) <- c("SLIDE", "SUBJECT", "ORGAN", "MI FINDING")

# Write to xlsx with each dataframe as a separate tab
write_xlsx(tab_list, "foresight_combined_study_data.xlsx")


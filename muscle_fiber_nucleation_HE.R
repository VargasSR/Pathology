setwd("~/Data Landing Zone/Muscle_Fiber")
getwd()
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

FILE <- "24IN0663_HE_Nucleation.tsv"
raw_table <- fread(FILE) 

#separate "name" by " " delimiter to make new columns
name_separated_table <- separate(raw_table, Name, c("study", "ID", "tissue", "group", "section", "date", "stain", "NA", "slide"), sep = " ")

ID_list <- unique(name_separated_table$ID)

#subset gastroc data
gastroc_table_names <- c("ID", "group", "section", grep("gastroc", names(name_separated_table), value=TRUE))

gastroc_nuclei_counts  <- na.omit(name_separated_table[, gastroc_table_names])

#subset soleus data
soleus_table_names <- c("ID", "group", "section", grep("Soleus", names(name_separated_table), value=TRUE))

soleus_nuclei_counts  <- na.omit(name_separated_table[, soleus_table_names])


# Gastroc Summary
gastroc_summary_table <- gastroc_nuclei_counts %>%
  group_by(ID, group, section) %>%
  summarize(
    pos_count = sum(`nuceli count gastroc` != 0),
    total_count = n(),
    percent_negative = 100 * (1 -(pos_count/ total_count))
  )

# soleus Summary
soleus_summary_table <- soleus_nuclei_counts %>%
  group_by(ID, group, section) %>%
  summarize(
    pos_count = sum(`nuceli count Soleus` != 0),
    total_count = n(),
    percent_negative = 100 * (1- (pos_count / total_count))
  )

# export
write.csv(soleus_summary_table, "soleus_summary_table.csv")
write.csv(gastroc_summary_table, "gastroc_summary_table.csv")
library(readxl)
library(writexl)
library(dplyr)

#Specify Directory
setwd("~/Data Landing Zone/Foresight_KD")
getwd()

df <- read_excel("Foresight_KD_20250718.xlsx", sheet = "MI FINDING") 

#remove Ds
df <- subset(df, KD != "D")

#remaining SUBJIDS
unique_subjids <- unique(df$SUBJID)

# Load SLIDE and SUBJECT tabs
df_slide <- read_excel("Foresight_KD_20250718.xlsx", sheet = "SLIDE")
df_subject <- read_excel("Foresight_KD_20250718.xlsx", sheet = "SUBJECT")

# Filter to only rows where SUBJID is in unique_subjids
df_slide <- df_slide[df_slide$SUBJID %in% unique_subjids,]
df_subject <- df_subject[df_subject$SUBJID %in% unique_subjids,]

# Create mapping from MISPEC to BlockID
  mapping <- c(
    "Kidney" = "A1",
    "Liver" = "B1", 
    "Lung" = "D1",
    "Trachea" = "D1",
    "Artery, Aorta" = "E1",
    "Heart" = "E1",
    "Gland, Parathyroid" = "F1",
    "Gland, Thyroid" = "F1",
    "Gland, Adrenal" = "F2",
    "Gland, Pituitary" = "F2",
    "Lymph Node, Mesenteric" = "G2",
    "Spleen" = "G2",
    "Thymus" = "G2"
  )

# Create new BlockID column
df$BLOCKID <- mapping[match(df$MISPEC, names(mapping))]

# Create VectorID column 
df$VectorID <- paste(df$SUBJID, df$BLOCKID, sep = "_")

# Create VectorID column
df_slide$VectorID <- paste(df_slide$SUBJID, df_slide$BLOCKID, sep = "_")

#unique vectorIDs
unique_vectorID <- unique(df$VectorID)

# Filter to only vectorIDs
df_slide <- df_slide[df_slide$VectorID %in% unique_vectorID,]

# Get duplicate VectorID rows 
#dup_rows <- df_slide[duplicated(df_slide$VectorID) | duplicated(df_slide$VectorID, fromLast = TRUE),]


# Create list of duplicate VectorIDs
#dup_vectorids <- unique(dup_rows$VectorID)

#save dupes
#write_xlsx(dup_rows, "Foresight_dupe_rows.xlsx")

# Delete columns
df$VectorID <- NULL 
df_slide$VectorID <- NULL


write_xlsx(df, "MI_FINDING.xlsx")
write_xlsx(df_slide, "SLIDE.xlsx") 
write_xlsx(df_subject, "SUBJECT.xlsx")
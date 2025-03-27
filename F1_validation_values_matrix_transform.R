setwd("~/Pathology/F1_validation/test_data/Averaged_F1_Data/3")
getwd()
library(data.table)
library(tidyr)

FILE <- "23IN0039 004 Ear A F 24Feb23 HE NA 4.tsv"
validation_raw_data <- fread(FILE) 

# Define column names
validation_values_column_names <- grep("ROI", names(validation_raw_data), value = TRUE)

# Subset to new dataframe with selected columns
validation_values <- na.omit(validation_raw_data[, ..validation_values_column_names])

# Extract the values from the single row
values <- as.numeric(validation_values[1, ])

# Reshape the values into a 10x10 matrix
matrix_values <- matrix(values, nrow = 10, byrow = TRUE)

# Assign column names (ROI1, ROI2, ..., ROI10)
colnames(matrix_values) <- paste0("ROI", 1:10)

# Assign row names (Label01, Label02, ..., Label10)
rownames(matrix_values) <- paste0("Label", sprintf("%02d", 1:10))

# Convert to data frame 
matrix_df <- as.data.frame(matrix_values)

# change names of ROIS and Labels 

colnames(matrix_df)[colnames(matrix_df) == "ROI1"] <- "Epidermis GT"
colnames(matrix_df)[colnames(matrix_df) == "ROI4"] <- "Muscle GT"
colnames(matrix_df)[colnames(matrix_df) == "ROI5"] <- "Dermis GT"
colnames(matrix_df)[colnames(matrix_df) == "ROI6"] <- "Adipose GT"
colnames(matrix_df)[colnames(matrix_df) == "ROI7"] <- "Cartilage GT"

rownames(matrix_df)[rownames(matrix_df) == "Label01"] <- "Epidermis Test"
rownames(matrix_df)[rownames(matrix_df) == "Label03"] <- "Muscle Test"
rownames(matrix_df)[rownames(matrix_df) == "Label04"] <- "Dermis Test"
rownames(matrix_df)[rownames(matrix_df) == "Label05"] <- "Adipose Test"
rownames(matrix_df)[rownames(matrix_df) == "Label06"] <- "Cartilage Test"

#  rename matrix to use in  published code

test01  <-  matrix_values

# define function
f1scores <- function(mat, conf.level=0.95){
  
  r <- ncol(mat)
  n <- sum(mat)
  
  p <- mat/n
  pii <- diag(p)
  pi. <- rowSums(p) 
  p.i <- colSums(p)
  
  miP <- miR <- sum(pii)
  miF1 <- miP
  
  F1i <- 2*pii/(pi.+p.i)
  maF1 <- sum(F1i)/r
  
  maP <- sum(pii/rowSums(p))/r
  maR <- sum(pii/colSums(p))/r
  maF2 <- 2*(maP*maR)/(maP+maR)
  
  miF1.v <- sum(pii)*(1-sum(pii))/n
  miF1.s <- sqrt(miF1.v)
  
  a <- 0
  b <- 0
  for(i in 1:r){
    jj <- (1:r)[-i]
    for(j in jj){
      b <- b + p[i,j]*F1i[i]*F1i[j]/((pi.[i]+p.i[i])*(pi.[j]+p.i[j]))
    }}
  maF1.v <- 2*(a+b)/(n*r^2) 
  maF1.s <- sqrt(maF1.v)
  
  varmap <- sum(pii*(pi.-pii)/pi.^3) / r^2 / n
  varmar <- sum(pii*(p.i-pii)/p.i^3) / r^2 / n
  covmpr1 <- sum( ((pi.-pii) * pii * (p.i-pii)) / (pi.^2 * p.i^2) )
  covmpr2 <- 0
  for(i in 1:r){
    covmpr2 <- covmpr2 + sum(pii[i] * p[i,-i] * pii[-i] / pi.[i]^2 / p.i[-i]^2)
  }
  covmpr <- (covmpr1+covmpr2) / r^2 / n
  maF2.v <- 4 * (maR^4*varmap + 2*maP^2*maR^2*covmpr + maP^4*varmar) / (maP+maR)^4
  maF2.s <- sqrt(maF2.v)
  
  z <- qnorm(1-(1-conf.level)/2)
  miF1.ci <- miF1 + c(-1,1)*z*miF1.s
  maF1.ci <- maF1 + c(-1,1)*z*maF1.s
  maF2.ci <- maF2 + c(-1,1)*z*maF2.s
  
  pr <- as.data.table(data.frame(microPrecision=miP, microRecall=miR, macroPrecision=maP, macroRecall=maR))
  
  fss <- as.data.table(data.frame(
    rbind(miF1=c(miF1, miF1.s, miF1.ci),  
          maF1=c(maF1, maF1.s, maF1.ci),
          maF1.star=c(maF2, maF2.s, maF2.ci))))
  
  setnames(fss, c('PointEst','Sd', 'Lower','Upper'))
  
  return(list(pr, fss))
}

pr <- f1scores(test01)[[1]]
fss <- f1scores(test01)[[2]]



# Create confusion matrix for visual
cm_table <- data.frame(test01) 

# Calculate the total for each ROI (columns)
total_ROI <- data.frame(t(apply(cm_table, 2, function(x) if(is.numeric(x)) sum(x) else "Total")))

# Calculate the total for each label (rows)
total_labels <- rowSums(cm_table[, sapply(cm_table, is.numeric)])
total_labels <-data.frame(total_labels)

# Convert total_ROI dataframe to a vector if necessary
total_ROI_vector <- as.vector(total_ROI)

# Create a new dataframe to store the results
cm_prop <- cm_table

# Loop through each element and perform the division with a check for zero
for (i in 1:nrow(cm_table)) {
  for (j in 1:ncol(cm_table)) {
    if (total_ROI_vector[j] != 0) {
      cm_prop[i, j] <- cm_table[i, j] / total_ROI_vector[j]
    } else {
      cm_prop[i, j] <- NA  
    }
  }
}

# Remove columns that contain only NA

cm_prop_clean <- cm_prop[, !apply(cm_prop, 2, function(x) all(is.na(x)))]

# remove rows that only contain zero
cm_prop_clean <- cm_prop_clean[!apply(cm_prop_clean, 1, function(x) all(x == 0)), ]

heatmap(as.matrix(cm_prop_clean), Rowv = NA, Colv = NA, scale = "none", col = hcl.colors(50))

cm_prop_log <- log10(cm_prop_clean + 1e-5)  # Adding a small constant to avoid log(0)

heatmap(as.matrix(cm_prop_log), Rowv = NA, Colv = NA, scale = "none", col = hcl.colors(50))



# Generate confusion matrix for calculation
cm <- as.table(test01)


# Initialize data frame 
class_results <- data.frame(Class = rownames(cm),
                            Precision = NA,
                            Recall = NA,
                            F1 = NA)

# Calculate metrics for each class
for(i in 1:nrow(cm)){
  tp <- cm[i,i]
  fp <- sum(cm[i,]) - tp
  fn <- sum(cm[,i]) - tp
  tn <- sum(cm) - tp - fp - fn
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * precision * recall / (precision + recall)
  
  class_results$Precision[i] <- precision
  class_results$Recall[i] <- recall
  class_results$F1[i] <- f1
}

# Calculate averaged metrics
avg_precision <- sum(diag(cm)) / sum(cm)
avg_recall <- sum(diag(cm)) / sum(rowSums(cm)) 
avg_f1 <- 2 * avg_precision * avg_recall / (avg_precision + avg_recall)

avg_results <- data.frame(Measure = c("Precision", "Recall", "F1"),
                          Value = c(avg_precision, avg_recall, avg_f1))

# Print tables  
print("Class-level metrics:")
print(class_results)

print("Average metrics:")
print(avg_results)

# Format results into table
results_table <- results
results_table$Measure <- rownames(results_table)
results_table <- reshape2::melt(results_table, id.vars="Measure")

# Print results table  
print("Results table:")
print(results_table)


# Calculate total annotations  
total_annotations <- sum(colSums(cm))

# Get proportion for each class
class_props <- colSums(cm) / total_annotations

# Calculate class-level metrics
class_precision <- diag(cm) / rowSums(cm)
class_recall <- diag(cm) / colSums(cm)
class_f1 <- 2 * class_precision * class_recall / (class_precision + class_recall)

# Calculate weighted averages
weighted_avg_precision <- sum(class_props * class_precision) 
weighted_avg_recall <- sum(class_props * class_recall)
weighted_avg_f1 <- sum(class_props * class_f1)





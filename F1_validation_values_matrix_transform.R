setwd("~/Data Landing Zone/Averaged F1 Validation")
getwd()
library(data.table)
library(tidyr)

FILE <- "test01_full_fill.tsv"
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

#  rename matrix for published code

mat <- matrix( matrix_values)


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
  
  a <- 0
  b <- 0
  for(i in 1:r){
    jj <- (1:r)[-i]
    for(j in jj){
      b <- b + p[i,j]*F1i[i]*F1i[j]/((pi.[i]+p.i[i])*(pi.[j]+p.i[j]))
    }}
  
  maF1.v <- 2*(a+b)/(n*r^2)
  maF1.s <- sqrt(maF1.v)
  
  miF1.v <- sum(pii)*(1-sum(pii))/n
  miF1.s <- sqrt(miF1.v)
  
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
  
  pr <- data.frame(microPrecision=miP, microRecall=miR, macroPrecision=maP, macroRecall=maR)
  fss <- data.frame(
    rbind(miF1=c(miF1, miF1.s, miF1.ci),
          maF1=c(maF1, maF1.s, maF1.ci),
          maF1.star=c(maF2, maF2.s, maF2.ci)))
  names(fss) <- c('PointEst','Sd', 'Lower','Upper')
  
  out <- list(pr, fss)
  names(out) <- c('Precision.and.Recall', 'Confidence.Interval')
  
  return(list(out, miF1.s, maF1.s))
  
  # Access standard errors
  mi_se <- results[[3]] 
  ma_se <- results[[4]]
  
}


results <- f1scores(mat)




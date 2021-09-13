

# directories
working_dir <- getwd()
setwd(working_dir)

# load libraries
library(ellipsenm)
library(kuenm)
library(usdm)
library(raster)
library(tidyverse)

# read layers
variables <- stack(list.files(path="sdm/calibration-areas/A_distichus-calibration-area", pattern=".asc", full.names=TRUE))
    plot(variables)
    # getting data from the variables
    variables_values <- na.omit(values(variables))

# sample of 50000 values if more pixels exist
if (dim(variables_values)[1] > 50000) {
  variables_values <- variables_values[sample(1:nrow(variables_values), 50000), ] 
}

# correlation matrix calculation
correlation_matrix <- cor(variables_values)

# saving correlation matrix
write.csv(correlation_matrix, "sdm/calibration-areas/A_distichus-calibration-area/variables_correlation_matrix.csv",
          row.names = TRUE)

# detecting correlated varaibles more easily
correlation_matrix1 <- correlation_matrix # making other table with results

max_cor <- 0.9 # maximum value of correlation allowed

for (i in 1:dim(correlation_matrix1)[2]) { #correlated values will turn into 2 for easier detection
  for (j in 1:dim(correlation_matrix1)[1]) {
    correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] > max_cor | correlation_matrix1[j, i] < -max_cor, 
                                        2, correlation_matrix1[j, i])
  }
}

# saving correlation matrix with rows exceeding maximum correlation threshold marked
write.csv(correlation_matrix1, "sdm/calibration-areas/A_distichus-calibration-area/variables_correlation_matrix_collinear_marked.csv",
          row.names = TRUE)

# #checking the table
View(correlation_matrix1) # selection should be done manually, 2 = correlated
names(variables)

usdm::vifcor(x=variables_values, th=0.85)

dir.create("sdm/non-correlated-variables")
selected_variables <- variables[[c()]]
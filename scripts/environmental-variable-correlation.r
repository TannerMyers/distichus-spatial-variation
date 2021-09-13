

# directories
working_dir <- getwd()
setwd(working_dir)

# load libraries
library(ellipsenm)
library(kuenm)
library(usdm)
library(corrplot)
library(raster)
library(tidyverse)

# read layers
variables <- raster::stack(list.files(path="sdm/calibration-areas/A_distichus-calibration-area", pattern=".asc", full.names=TRUE))
    # plot(variables)
    # getting data from the variables
    variables_values <- na.omit(values(variables))

# sample of 50000 values if more pixels exist
if (dim(variables_values)[1] > 50000) {
  variables_values <- variables_values[sample(1:nrow(variables_values), 50000), ] 
}

# correlation matrix calculation
correlation_matrix <- cor(variables_values)

  ## visualize relationships between variables
  corrplot(correlation_matrix)
  corrplot.mixed(correlation_matrix, upper="ellipse", lower="number")

# saving correlation matrix
write.csv(correlation_matrix, "sdm/calibration-areas/A_distichus-calibration-area/variables_correlation_matrix.csv",
          row.names = TRUE)

# detecting correlated varaibles more easily
correlation_matrix1 <- correlation_matrix # making other table with results

max_cor <- 0.85 # maximum value of correlation allowed

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


# Estimate Variance Influence Factors (VIFs) for each of the predictors
  ## `vifcor` function of the `usdm` package estimates VIFs for a Raster Stack
  ## by finding a pair of variables with a maximum linear correlation and excludes
  ## the one with the greater VIF until no variable with a high correlation coefficient
  ## with other variables remains
usdm::vifcor(x=variables_values, th=0.85) # 0.85 is linear correlation

  ## `vifstep` function of the `usdm` package estimates VIFs for a Raster Stack
  ## all at once, excluding the one with the highest VIF, repeating until no variables
  ## with VIF higher than th remains
usdm::vifstep(x=variables_values, th=10) # 10 is the threshold value of VIF 

# Create a directory to contain non-correlated variables
dir.create("sdm/non-correlated-variables")

selected<- c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_05", 
                                  "CHELSA_bio10_15", "CHELSA_bio10_16", "march_EVI_mean",
                                  "may_NDVI_mean")

selected_variables <- variables[[selected]]

variable_names <- paste0("sdm/non-correlated-variables/", selected, ".asc")

for (i in 1:nlayers(selected_variables)) {
  writeRaster(selected_variables[[i]], filename= variable_names[i], format="ascii",)
}

# Now, output combinations of selected variables for Maxent ENM estimation with `kuenm` package
vs <- kuenm_varcomb(var.dir = "sdm/non-correlated-variables", out.dir="sdm/M_variables", min.number=3, in.format ="ascii", out.format="ascii")
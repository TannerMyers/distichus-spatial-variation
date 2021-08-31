
# directories
working_dir <- getwd()
setwd(working_dir)

# load libraries
library(ellipsenm)
library(raster)

# read layers
variables <- stack(list.files())
working_dir <- getwd()
setwd(working_dir)

# Load R packages
library(ENMeval)
library(raster)
library(dplyr)

set.seed()

occ_dir <- paste0(working_dir, "/data/occurrences/")
occs <- read_csv(paste0(occ_dir, "filtered_distichus_occurrences.csv"))
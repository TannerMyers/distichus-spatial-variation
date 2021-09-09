
# directories
working_dir <- getwd()
setwd(working_dir)

# load libraries
library(ellipsenm)
library(raster)
library(tidyverse)

# read layers
variables <- stack(list.files(path="sdm/calibration-areas/A_distichus-calibration-area", pattern=".asc", full.names=TRUE))

# load records 
records <- read_csv("sdm/thinned-datasets/A_distichus_30km.csv")
    ## consists of subspecies names as is 
    records$scientificName <- "Anolis distichus"

# Exploring environmental space --------------------------------------------------------

# exploring the data in environmental space
## only temperature variables
explore_espace(data = records, species ="scientificName", longitude = "longitude",
                latitude = "latitude", raster_layers = variables[[c(1:9,20)]], save = TRUE, 
                name = "temp_plots/Temperature_Elevation_variables.pdf")

## only precipitation variables
explore_espace(data = records, species ="scientificName", longitude = "longitude",
                latitude = "latitude", raster_layers = variables[[10:19]], save = TRUE, 
                name = "temp_plots/Precipitation_VI_variables.pdf")

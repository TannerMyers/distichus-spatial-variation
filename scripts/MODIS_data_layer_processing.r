###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# This script 

###############

# load packages
library(raster)
library(rgdal)
library(rgeos)
library(tidyverse)

# Load rasters for May vegetation indices
    ## Files should have '121' following the year. This is the Julian date for May 1st
may_evi <- raster::stack(list.files(path="May", pattern="EVI", full.names = TRUE, recursive = TRUE))
may_ndvi <- raster::stack(list.files(path="May", pattern="NDVI", full.names = TRUE, recursive = TRUE))

# Load rasters for March vegetation indices
    ## Files should have '060' following the year. This is the Julian date for March 1st
march_evi <- raster::stack(list.files(path="March", pattern="EVI", full.names = TRUE, recursive = TRUE))
march_ndvi <- raster::stack(list.files(path="March", pattern="NDVI", full.names = TRUE, recursive = TRUE))



may_evi_mean <- raster::calc(may_evi, mean)
may_ndvi_mean <- raster::calc(may_ndvi, mean)

march_evi_mean <- raster::calc(march_evi, mean)
march_ndvi_mean <- raster::calc(march_ndvi, mean)

# Save mean rasters to files
writeRaster(x = march_ndvi_mean, filename =file.path("March/march_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_ndvi_mean, filename = file.path("May/may_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = march_evi_mean, filename = file.path("March/march_EVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_evi_mean, filename = file.path("May/may_EVI_mean"), format = "ascii", overwrite = TRUE)


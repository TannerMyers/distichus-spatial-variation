###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# This script averages yearly data layers for the MODIS variables NDVI and
# EVI for the months of March and May to create mean monthly EVI and NDVI
# for the decade of interest.

###############

# directories 
working_dir <- getwd()
setwd(working_dir)

# load packages
library(raster)
library(rgdal)
library(rgeos)
library(tidyverse)

# Load rasters for May vegetation indices
    ## Files should have '121' following the year. This is the Julian date for May 1st
may_evi <- raster::stack(list.files(path="data/MODIS/May", pattern="monthly_EVI", full.names = TRUE, recursive = TRUE))
may_ndvi <- raster::stack(list.files(path="data/MODIS/May", pattern="monthly_NDVI", full.names = TRUE, recursive = TRUE))

# Load rasters for March vegetation indices
    ## Files should have '060' following the year. This is the Julian date for March 1st
march_evi <- raster::stack(list.files(path="data/MODIS/March", pattern="monthly_EVI", full.names = TRUE, recursive = TRUE))
march_ndvi <- raster::stack(list.files(path="data/MODIS/March", pattern="monthly_NDVI", full.names = TRUE, recursive = TRUE))


may_evi_mean <- raster::calc(may_evi, mean)
may_ndvi_mean <- raster::calc(may_ndvi, mean)

march_evi_mean <- raster::calc(march_evi, mean)
march_ndvi_mean <- raster::calc(march_ndvi, mean)

dir.create("data/MODIS_new")

# Save mean rasters to files
writeRaster(x = march_ndvi_mean, filename =file.path("data/MODIS_new/march_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_ndvi_mean, filename = file.path("data/MODIS_new/may_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = march_evi_mean, filename = file.path("data/MODIS_new/march_EVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_evi_mean, filename = file.path("data/MODIS_new/may_EVI_mean"), format = "ascii", overwrite = TRUE)


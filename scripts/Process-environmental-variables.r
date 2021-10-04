###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# Code modified from https://github.com/townpeterson/vespa/tree/master/Rcode

# This script loads and processes bioclimatic variables from the CHELSA database,
# which I downloaded previously using wget on the urls present in "data/envidatS3paths.txt"

##############

# directories 
working_dir <- getwd()
setwd(working_dir)

# libraries
library(raster)
library(rgdal)
library(rgeos)
library(tidyverse)
library(httpgd)
#----------------------------------------------------------------

# Process bioclimatic variables from CHELSA database ------------
chelsa <- stack(list.files("data/chelsa", pattern=".tif", full.name=TRUE))
hgd()
hgd_browse()
plot(chelsa)

        ## I am not interested in projecting the suitability of habitats outside of Hispaniola 
        ## so I am going to crop both the raster & shapefile loaded below
        limits <- extent(c(-76,-67, 17, 21))
        chelsa <- crop(chelsa, limits)

# Shapefile featuring country boundaries downloaded from: https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0/data
    world <- readOGR(dsn = ".", "World_Countries__Generalized_")
    world$id <- "1"
    hisp <- crop(world, limits)

    nhisp <- gUnaryUnion(hisp, hisp@data$id)

    # Getting the cells with values 
    cells_nona <- cellFromPolygon(chelsa[[1]], nhisp )
    cellsids <- 1:ncell(chelsa[[1]])
    cell_na <- cellsids[-cells_nona[[1]]]
    chelsa[cell_na] <- NA
    plot(chelsa[[1]])

# Saving variables to a new directory
    dir.create("data/chelsa_new")

    lapply(names(chelsa), function(x) {
        writeRaster(chelsa[[x]], paste0("data/chelsa_new/", x, ".asc"), overwrite = TRUE)
    })

#----------------------------------------------------------------

# Process MODIS variables ---------------------------------------

## Load rasters for May vegetation indices
    ## Files should have '121' following the year. This is the Julian date for May 1st
may_evi <- raster::stack(list.files(path="data/MODIS/May", pattern="monthly_EVI", full.names = TRUE, recursive = TRUE))
may_ndvi <- raster::stack(list.files(path="data/MODIS/May", pattern="monthly_NDVI", full.names = TRUE, recursive = TRUE))

## Load rasters for March vegetation indices
    ## Files should have '060' following the year. This is the Julian date for March 1st
march_evi <- raster::stack(list.files(path="data/MODIS/March", pattern="monthly_EVI", full.names = TRUE, recursive = TRUE))
march_ndvi <- raster::stack(list.files(path="data/MODIS/March", pattern="monthly_NDVI", full.names = TRUE, recursive = TRUE))

## estimate mean for each layer that represent a year
may_evi_mean <- raster::calc(may_evi, mean)
may_ndvi_mean <- raster::calc(may_ndvi, mean)

march_evi_mean <- raster::calc(march_evi, mean)
march_ndvi_mean <- raster::calc(march_ndvi, mean)

## create directory to output averaged rasters to files
dir.create("data/MODIS_new")

## Save mean rasters to files
writeRaster(x = march_ndvi_mean, filename =file.path("data/MODIS_new/march_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_ndvi_mean, filename = file.path("data/MODIS_new/may_NDVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = march_evi_mean, filename = file.path("data/MODIS_new/march_EVI_mean"), format = "ascii", overwrite = TRUE)
writeRaster(x = may_evi_mean, filename = file.path("data/MODIS_new/may_EVI_mean"), format = "ascii", overwrite = TRUE)

#----------------------------------------------------------------

# Process elevation data from SRTM database ---------------------

## Include elevation among other environmental variables
dir.create("data/elevation")
raster::getData(name="SRTM", path="data/elevation", lon = -74, lat=18)
raster::getData(name="SRTM", path="data/elevation", lon = -69, lat=18)
raster::getData(name="SRTM", path="data/elevation", lon = -74, lat=21)

## load downloaded SRTM tiles
alt1 <- raster("data/elevation/srtm_22_09.tif")
alt2 <- raster("data/elevation/srtm_23_09.tif")
alt3 <- raster("data/elevation/srtm_22_08.tif")

## Need to merge rasters because they represent different tiles from the SRTM dataset
elev <- merge(alt1,alt2,alt3)

## crop to match MODIS extent, using the raster for mean March EVI
elev <- crop(elev, march_evi_mean@extent)
    elev@extent <- raster::alignExtent(elev@extent, march_evi_mean) 

    ## we need to resample the elevation raster because it has 
    ## 10x as many rows and columns than CHELSA and MODIS variables
    ## because SRTM is at ~90 m resolution, while CHELSA and MODIS are ~1 km
    elev_resampled <- resample(elev, march_evi_mean)

    ## write directory to output resampled elevation raster
    dir.create("data/elevation_new")
    ## Now, output the resampled raster
    writeRaster(x = elev_resampled, filename =file.path("data/elevation_new/SRTM_elevation_1km"), format = "ascii", overwrite = TRUE)
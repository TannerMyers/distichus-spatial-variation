
# directories 
working_dir <- getwd()
setwd(working_dir)

# libraries
library(raster)
library(rgeos)
library(rgdal)
library(RStoolbox)
library(tidyverse)


bc_list <- raster::stack(list.files("data/chelsa_new/", pattern=".asc", full.names = TRUE))
vi_list <- raster::stack(list.files("data/MODIS_new/", pattern=".asc", full.names = TRUE))

## Include elevation among other environmental variables
    # dir.create("data/elevation")
    # raster::getData(name="SRTM", path="data/elevation", lon = -74, lat=18)
    # raster::getData(name="SRTM", path="data/elevation", lon = -69, lat=18)
    # raster::getData(name="SRTM", path="data/elevation", lon = -74, lat=21)
alt1 <- raster("data/elevation/srtm_22_09.tif")
alt2 <- raster("data/elevation/srtm_23_09.tif")
alt3 <- raster("data/elevation/srtm_22_08.tif")
## Need to merge rasters because they represent different tiles from the SRTM
## dataset
elev <- merge(alt,alt2,alt3)

    ## Rasters have different extents
        extent(bc_list)
        # class      : Extent 
        # xmin       : -76.00014 
        # xmax       : -67.00014 
        # ymin       : 16.99986 
        # ymax       : 20.99986 
        extent(vi_list)
        # class      : Extent 
        # xmin       : -74.81667 
        # xmax       : -68.16667 
        # ymin       : 17.4 
        # ymax       : 20.81667 
    ## Since the extent of vi_list fits within bc_list's extent, let's
    ## crop bc_list to make it the same for both
    bc_list <- crop(bc_list, vi_list@extent)
    ## still not quite the same
    bc_list@extent <- raster::alignExtent(bc_list@extent, vi_list)

    elev <- crop(elev, vi_list@extent)
    elev@extent <- raster::alignExtent(elev@extent, vi_list) 

    ## we need to resample the elevation raster because it has 
    ## 10x as many rows and columns than CHELSA and MODIS variables
        ## SRTM is at ~90 m resolution, while CHELSA and MODIS are ~1 km
    elev_resampled <- resample(elev, vi_list)


## Now, merge resolved raster stacks 
env <- raster::stack(elev_resampled, bc_list, vi_list, full.names=TRUE)
plot(env)

# 

## perform principal component analysis on rasters
pca <- rasterPCA(env, spca=TRUE)
    ## Visualize PC rasters
    plot(pca$map)

## Check loadings & eigenvalues for top 5 loadings
knitr::kable(round(pca$model$loadings[,1:5],5)) # loadings
summary(pca$model) # eigenvalues

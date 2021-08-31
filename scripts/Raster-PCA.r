
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
elev <- raster("data/elevation_new/SRTM_elevation_1km.asc")

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

## Now, merge resolved raster stacks 
env <- raster::stack(elev, bc_list, vi_list, full.names=TRUE)
plot(env)

## perform principal component analysis on rasters
pca <- rasterPCA(env, spca=TRUE)
    ## Visualize PC rasters
    plot(pca$map)

## Check loadings & eigenvalues for top 5 loadings
knitr::kable(round(pca$model$loadings[,1:5],5)) # loadings
summary(pca$model) # eigenvalues

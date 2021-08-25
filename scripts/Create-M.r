
# directory
working_dir <- getwd()
setwd(working_dir)

# load libraries 
library(raster)
library(rgeos)
library(rgdal)
library(maps)

# Setting variables

## considering earth's distortion. 
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# define function to create M 
create_M <- function(occ) {
    occ_sp <- SpatialPointsDataFrame(coords=occ[,2:3], data=occ, proj4string=WGS84)

    ## project points using their centroids as a reference
}


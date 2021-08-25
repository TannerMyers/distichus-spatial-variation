


# load libraries 
library(raster)
library(rgeos)
library(rgdal)
library(maps)
library(tidyverse)

# Define function create_M
    ## Takes two inputs: a dataframe with cleaned records for your taxon of choice and
    ## the taxon's name 
create_M <- function(occ, taxon) {
    occ_sp <- SpatialPointsDataFrame(coords=occ[,2:3], data=occ, proj4string=WGS84)

    ## project points using their centroids as a reference
    centroid <- gCentroid(occ_sp, byid = FALSE)
    AEQD <- CRS(paste("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", centroid@coords[1],
                  " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", sep = ""))
    
    occ_proj <- spTransform(occ_sp, AEQD)

    ## create a buffer based on 50 km distance 
    buff_area <- gBuffer(occ_proj, width = 50000, quadsegs = 30)
    buff_area <- disaggregate(buff_area)

    ## reproject buffered area
    buff_area <- spTransform(buff_area, WGS84)

    ## make spatialpolygondataframe
    df <- data.frame(species = rep(taxon, length(buff_area)))
    buff_area <- SpatialPolygonsDataFrame(buff_area, data = df, match.ID = FALSE)

    ## write area out as a shapefile
    shp_dir <- paste0(output_dir,taxon,"-calibration-area")
    dir.create(shp_dir)
    writeOGR(buff_area, shp_dir, "M", driver = "ESRI Shapefile")

    ## plot 
    lims <- extent(buff_area)
    map(xlim = lims[1:2] + c(-1,1), ylim = lims[3:4] + c(-1,1))
    points(buff_area, border = "purple", add = TRUE, lwd = 2)
}

# Setting variables

## directory
working_dir <- getwd()
setwd(working_dir)

# create directory to output files
output_dir <- "enm/"
dir.create(output_dir)

## considering earth's distortion. 
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

## Load inputs for `create_M`
all_distichus <- read_csv("data/GBIF/A_distichus_clean.csv")
distichus <- "A_distichus"

# 
create_M(all_distichus, distichus)

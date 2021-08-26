###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# Code modified from https://github.com/townpeterson/vespa/tree/master/Rcode

# This script creates "M" for ecological suitability modeling and masks bioclimatic variables
# from the calibration areas it creates.


# load libraries 
library(raster)
library(rgeos)
library(rgdal)
library(maps)
library(tidyverse)

## directory
working_dir <- getwd()
setwd(working_dir)

#---------------------------------------------------------------------

#---------------- Preparing calibration areas ------------------------


# Setting variables

# create directory to output files
output_dir <- "sdm/calibration-areas/"
dir.create(output_dir)

## considering earth's distortion. 
WGS84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Define function `create_M`
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
    maps::map(xlim = lims[1:2] + c(-1,1), ylim = lims[3:4] + c(-1,1))
   # points(buff_area, border = "purple", add = TRUE, lwd = 2)
   ## ^Resulted in following error:
   ## Error in as.double(y) : 
   ## cannot coerce type 'S4' to vector of type 'double'

}

#

## Load inputs for `create_M`
all_distichus_occs <- read_csv("sdm/thinned-datasets/A_distichus_30km.csv")
distichus <- "A_distichus"
dom_occs <- read_csv("sdm/thinned-datasets/A_d_dominicensis_30km.csv")
dom <- "A_d_dominicensis"
ign_occs <- read_csv("sdm/thinned-datasets/A_d_ignigularis_10km.csv")
ign <- "A_d_ignigularis"
rav_occs <- read_csv("sdm/thinned-datasets/A_d_ravitergum_10km.csv")
rav <- "A_d_ravitergum"
fav_occs <- read_csv("sdm/thinned-datasets/A_d_favillarum_3km.csv")
fav <- "A_d_favillarum"
prop_occs <- read_csv("sdm/thinned-datasets/A_d_properus_10km.csv")
prop <- "A_d_properus"
haiti_occs <- read_csv("sdm/thinned-datasets/A_distichus_Tiburon_10km.csv")
haiti <- "A_distichus_Tiburon"

## Run it on Anolis distichus complete dataset
create_M(all_distichus_occs, distichus)
create_M(dom_occs, dom)
create_M(fav_occs, fav)
create_M(haiti_occs, haiti)
create_M(ign_occs, ign)
create_M(prop_occs, prop)
create_M(rav_occs, rav)

#-------------------------------------------------------------------------------


#---------------------Masking variables to calibration area---------------------

## load bioclimatic variables
chelsa_clim <- raster::stack(list.files(path = "data/chelsa_new/", pattern = ".asc", full.names = TRUE))

dirs <- list.files("sdm/calibration-areas/")
for (dir in dirs){
    dir <- str_remove(dir, '[1]')
    M <- readOGR(dir, layer = "M")

    varsm <- mask(crop(chelsa_clim, M), M)

    lapply(names(varsm), function(x) {
        writeRaster(varsm[[x]], paste0(dir,"/",x,".asc"), overwrite = TRUE)
    })
}

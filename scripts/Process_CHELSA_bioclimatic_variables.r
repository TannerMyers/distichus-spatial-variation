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

library(raster)
library(rgdal)
library(rgeos)
library(httpgd)

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

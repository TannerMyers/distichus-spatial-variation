###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# Code obtained and modified from https://github.com/townpeterson/vespa/tree/master/Rcode

# This script downloads raw MERRAclim data from Dryad and processes the bioclimatic variables
# for use in niche modeling with kuenm 

###############

# load packages
    library(raster)
    library(rgdal)
    library(rgeos)

#----------------------------Downloading & unzipping bioclimatic variables----------------------------#
        ## https://datadryad.org/stash/dataset/doi:10.5061/dryad.s2v81 is the MERRAclim data link
    ## download url for MerraClim
    data_url <- "https://datadryad.org/stash/downloads/file_stream/97913" # 2000s means @ 2.5 arc-minutes
    
    ## download data
    download.file(data_url, destfile = "2_5m_mean_00s.zip", method = "auto", 
              quiet = FALSE, mode = "wb", cacheOK = TRUE)

    ## unzip variables
    var_dir <- "2_5m_mean_00s"
    dir.create(var_dir)
    unzip(zipfile = "2_5m_mean_00s.zip", exdir=var_dir)
#-----------------------------------------------------------------------------------------------------#

#--------------------------------------Preparing variables--------------------------------------------#
    ## stacking of variables (excluding 8,9,18, and 19; see Escobar et al. (2014))
    mc2_5 <- stack(list.files("2_5m_mean_00s", pattern = "bio.*if$", full.names = T))
    names(mc2_5) <- gsub("X2_5m_mean_00s_", "", names(mc2_5))
    mc2_5 <- mc2_5[[c(1, 12:17, 2:9)]]

        ## At 2.5 arc-minutes, this is a lot of data to deal with and I am not interested in projecting
        ## the suitability of habitats outside of Hispaniola so I am going to crop both the raster & shapefile loaded below
        limits <- extent(c(-76,-67, 17, 21))
        mc2_5_subset <- crop(mc2_5, limits)

    # Shapefile featuring country boundaries from: https://hub.arcgis.com/datasets/252471276c9941729543be8789e06e12_0/data
    world <- readOGR(dsn = ".", "World_Countries__Generalized_")
    world$id <- "1"
    hisp <- crop(world, limits)

    nhisp <- gUnaryUnion(hisp, hisp@data$id)

    # Getting the cells with values 
    cells_nona <- cellFromPolygon(mc2_5_subset[[1]], nhisp )
    cellsids <- 1:ncell(mc2_5_subset[[1]])
    cell_na <- cellsids[-cells_nona[[1]]]
    mc2_5_subset[cell_na] <- NA
    plot(mc2_5_subset[[1]])

    # Saving variables to a new directory
    dir.create("merra_new")

    lapply(names(mc2_5_subset), function(x) {
        writeRaster(mc2_5_subset[[x]], paste0("merra_new/", x, ".asc"), overwrite = TRUE)
    })

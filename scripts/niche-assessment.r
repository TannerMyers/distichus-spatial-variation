rm(list = ls())

# Set working directory
working_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation"
setwd(working_dir)

# Load libraries
library(spocc)
library(sp)
library(raster)
library(rgeos)
library(sf)
library(rmapshaper)
library(envirem)
library(rangeBuilder)
library(fauxcurrence)
library(ENMTools)
library(usdm)
library(RStoolbox)
library(data.table)
library(tidyverse)

# Directory with filtered and cleaned occurrences
occ_dir <- paste0(working_dir, "/data/occurrences/")
    ## Load filtered occurrences downloaded with `spocc` R package
    occs_ad <- read_csv(paste0(occ_dir, "filtered_distichus_occurrences.csv"))
    occs_ab <- read_csv(paste0(occ_dir, "filtered_brevirostris_occurrences.csv"))

# Omit samples not from Hispaniola
    ## To do this, we will load a shapefile of the island of Hispaniola
    DOM <- getData('GADM', country = 'DOM', level=0, path = paste0(working_dir, "/data/shape-files/"), download = FALSE)
    HTI <- getData('GADM', country = 'HTI', level=0, path = paste0(working_dir, "/data/shape-files/"), download = FALSE)
    row.names(DOM) <- paste("DOM", row.names(DOM), sep = "_")
    row.names(HTI) <- paste("HTI", row.names(HTI), sep = "_")
    Hispaniola <- rbind(HTI, DOM, makeUniqueIDs = TRUE)
    Hispaniola <- gSimplify(Hispaniola, tol = 0.01, topologyPreserve = TRUE)

    ## Get extent for cropping
    limits <- extent(Hispaniola)


# Load environmental data
chelsa_clim <- raster::stack(list.files(path = "data/chelsa_new/", pattern = ".asc", full.names = TRUE))
modis_vi <- raster::stack(list.files(path = "data/MODIS_new/", pattern = ".asc", full.names = TRUE))
## These raster stacks have different extents
    ## Since the extent of modis_vi fits within chelsa_clim's extent, crop chelsa_clim
    chelsa_clim <- crop(chelsa_clim, modis_vi@extent)
    chelsa_clim@extent <- raster::alignExtent(chelsa_clim@extent, modis_vi)

## Now, load elevational data that has already been adjusted to the MODIS variable extent
elev_srtm <- raster("data/elevation_new/SRTM_elevation_1km.asc")

## Merge all environmental data layers into single raster stack
env <- raster::stack(elev_srtm, chelsa_clim, modis_vi, full.names = TRUE)
    env <- crop(env, Hispaniola)
    ## Check rasters for same extent and resolution &
    ## set values to NA if NAs present in any layer
    env <- ENMTools::check.env(env)

env2 <- env[[c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_05",
            "CHELSA_bio10_15", "CHELSA_bio10_16",
            "march_EVI_mean", "may_NDVI_mean")]]

# Load ddRAD sampling spreadsheet
rad_data <- read_table("~/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
    #rad_data[36,2] <- "1352_dom2" ## Note there is a discrepancy in how this specimen is classified between my sampling spreadsheets
cluster <- read_table("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/lea/distichus-complex-only/sNMF_K5_ancestry_coefficients_cluster-specimen-cleaned.tsv", col_names = TRUE)

rad_data <- rad_data[rad_data$Sample_ID_pop %in% cluster$Sample_ID_pop, ]
    rad_data <- as_tibble(cbind(rad_data, cluster$cluster))
    colnames(rad_data)[12] <- "snmfCluster" 
    ## Remove duplicates
    rad_data <- as_tibble(unique(setDT(rad_data), by = 'Locality'))

# Make range polygons and raster and sort occurrences by lineage
all_pts <- SpatialPoints(occs_ad[, 2:3])
    crs(all_pts) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

for (k in 1:5){
    pts <- SpatialPoints(rad_data[rad_data$snmfCluster == k, c("Longitude", "Latitude")])
        crs(pts) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

    poly <- rangeBuilder::getDynamicAlphaHull(x = pts@coords, fraction = 1.0,
                    clipToCoast = "terrestrial",
                    buff = 1000)
        assign(paste0("K", k, "_alpha_poly"), poly)
    
    lims <- extent(poly[[1]])
        assign(paste0("K", k, "_extent"), lims)

    ## Subset occurrences that fall within spatial polygon made above
    df <- occs_ad[!is.na(over(all_pts, poly[[1]])), ]
    # Combine total points
    colnames(pts@coords) <- colnames(all_pts@coords)
    df2 <- as_tibble(rbind(pts@coords, df[, 2:3]))
    df3 <- as_tibble(base::as.data.frame(cbind(df2, rep(k, nrow(df2)))))
        colnames(df3)[3] <- "K"
        df3$K <- as.factor(df3$K)
    write_delim(x = df3, file = paste0("niche-assessment/", "K", k, "_total_ah_pts.csv"),
                delim = ",", col_names = TRUE)
        assign(paste0("K", k, "_total_pts"), df3)
    
    ## Create raster for background point sampling
    crop <- crop(env2[[1]], lims)
    mask <- mask(crop, poly[[1]])
    crs(mask) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        writeRaster(x = mask, filename = paste0("niche-assessment/", "K", k, "_masked.asc"), format = "ascii", overwrite = TRUE)
        assign(paste0("K", k, "_bias"), mask)
}

## The second population has occurrences on satellite islands so let's prune those
K2_alpha_poly[[1]] <- ms_filter_islands(input = K2_alpha_poly[[1]], min_area = 1500000000)
    K2_spat <- SpatialPoints(K2_total_pts[, 1:2])
    crs(K2_spat) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
K2_total_pts <- K2_total_pts[!is.na(over(K2_spat, K2_alpha_poly[[1]])), ]
    write_delim(x = K2_total_pts, file = "niche-assessment/K2_total_ah_pts.csv", delim = ",", col_names = TRUE)
K2_extent <- extent(K2_alpha_poly[[1]])
crop <- crop(env2[[1]], K2_extent)
K2_bias <- mask(crop, K2_alpha_poly[[1]])
    crs(K2_bias) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    writeRaster(x = K2_bias, filename = "niche-assessment/K2_masked.asc", format = "ascii", overwrite = TRUE)

# Merge into one big dataframe
all_pts <- as_tibble(rbind(K1_total_pts, K2_total_pts, K3_total_pts, K4_total_pts, K5_total_pts))
    bg_input <- all_pts[, 1:3]
    colnames(bg_input) <- c("x", "y", "species")

# Use fauxcurrence to obtain bias-corrected background points
set.seed(69)
    ## Within-species distances i.e., "intra" model
    faux_intra <- fauxcurrence(coords = as.data.frame(bg_input), rast = env2[[1]])
        save(faux_intra, file = "niche-assessment/faux_intra.RData")
    ## Pairwise between species
    faux_interSep <- fauxcurrence(coords = as.data.frame(bg_input), rast = env2[[1]],
                            inter.spp = TRUE, sep.inter.spp = TRUE, verbose = TRUE)
        save(faux_interSep, file = "niche-assessment/faux_interSep.RData")
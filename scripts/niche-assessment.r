rm(list = ls())

# Set working directory
working_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation"
setwd(working_dir)

# Load libraries
library(sp)
library(raster)
library(rgeos)
library(sf)
library(rmapshaper)
library(envirem)
library(rangeBuilder)
library(ellipsenm)
library(fauxcurrence)
library(ENMTools)
library(usdm)
library(RStoolbox)
library(data.table)
library(tidyverse)

# Add imputed rasters to a new raster stack
env_imputed <- raster::stack(list.files(path = "data/env_imputed", pattern = ".asc", full.names = TRUE))

# Subset rasters left after tests for collinearity
env2 <- env_imputed[[c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_05",
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

    ## Create raster for background point sampling
    crop <- crop(env2[[1]], lims)
    mask <- mask(crop, poly[[1]])
    crs(mask) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        writeRaster(x = mask, filename = paste0("niche-assessment/", "K", k, "_masked.asc"), format = "ascii", overwrite = TRUE)
        assign(paste0("K", k, "_bias"), mask)

    ## Subset occurrences that fall within spatial polygon made above
    df <- occs_ad[!is.na(over(all_pts, poly[[1]])), ]
    # Combine total points
    colnames(pts@coords) <- colnames(all_pts@coords)
    df2 <- as_tibble(rbind(pts@coords, df[, 2:3]))
    df3 <- as_tibble(base::as.data.frame(cbind(df2, rep(k, nrow(df2)))))
        colnames(df3)[3] <- "K"
        df3$K <- as.factor(df3$K)
        df3 <- df3[!duplicated(paste(df3$longitude, df3$latitude)), ]
    write_delim(x = df3, file = paste0("niche-assessment/", "K", k, "_total_ah_pts.csv"),
                delim = ",", col_names = TRUE)
        assign(paste0("K", k, "_total_pts"), df3)
    
    ## Thin points
    df4 <- thin_data(df3, "longitude", "latitude", thin_distance = 15, save = FALSE)
        write_delim(x = df4, file = paste0("niche-assessment/", "K", k, "thinned.csv"), delim = ",", col_names = TRUE)
        assign(paste0("K", k, "_thinned"), df4)
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

# Prepare data for fauxcurrence    
colnames(all_pts) <- c("x", "y", "species")
all_pts$species <- as.character(all_pts$species)
all_pts$species <- replace(all_pts$species, all_pts$species == 1, "dominicensis")
all_pts$species <- replace(all_pts$species, all_pts$species == 2, "South")
all_pts$species <- replace(all_pts$species, all_pts$species == 3, "ignigularis")
all_pts$species <- replace(all_pts$species, all_pts$species == 4, "ravitergum")
all_pts$species <- replace(all_pts$species, all_pts$species == 5, "properus")

K_12 <- merge(K1_bias, K2_bias)
K_13 <- merge(K1_bias, K3_bias)
K_24 <- merge(K2_bias, K4_bias)
K_34 <- merge(K3_bias, K4_bias)
K_35 <- merge(K3_bias, K5_bias)

# Use fauxcurrence to obtain bias-corrected background points
set.seed(69)
    # # Pairwise between species pairs I will be performing niche similarity tests on
    # coords <- as.data.frame(all_pts[all_pts$species == "dominicensis" | all_pts$species == "South", ])
    # faux_inter_12 <-  fauxcurrence(coords = coords, rast = K_12, inter.spp = TRUE, sep.inter.spp = TRUE, allow.ident.conspec = FALSE)
    #     save(faux_inter_12, file = "niche-assessment/faux_inter_12.RData")

    # coords <- as.data.frame(all_pts[all_pts$species == "dominicensis" | all_pts$species == "ignigularis", ])
    # faux_inter_13 <-  fauxcurrence(coords = coords, rast = K_13, inter.spp = TRUE, sep.inter.spp = TRUE, allow.ident.conspec = FALSE)
    #     save(faux_inter_13, file = "niche-assessment/faux_inter_13.RData")

    # coords <- as.data.frame(all_pts[all_pts$species == "South" | all_pts$species == "ravitergum", ])
    # faux_inter_24 <- fauxcurrence(coords = coords, rast = K_24, inter.spp = TRUE, sep.inter.spp = TRUE, allow.ident.conspec = FALSE)
    #     save(faux_inter_24, file = "niche-assessment/faux_inter_24.RData")

    # coords <- as.data.frame(all_pts[all_pts$species == "ignigularis" | all_pts$species == "ravitergum", ])
    # faux_inter_34 <- fauxcurrence(coords = coords, rast = K_34, inter.spp = TRUE, sep.inter.spp = TRUE, allow.ident.conspec = FALSE)
    #     save(faux_inter_34, file = "niche-assessment/faux_inter_34.RData")

    # coords <- as.data.frame(all_pts[all_pts$species == "ignigularis" | all_pts$species == "properus", ])
    # faux_inter_35 <- fauxcurrence(coords = coords, rast = K_35, inter.spp = TRUE, sep.inter.spp = TRUE, allow.ident.conspec = FALSE)
    #     save(faux_inter_35, file = "niche-assessment/faux_inter_35.RData")
    
load("niche-assessment/faux_inter_12.RData")
load("niche-assessment/faux_inter_13.RData")
load("niche-assessment/faux_inter_24.RData")
load("niche-assessment/faux_inter_34.RData")
load("niche-assessment/faux_inter_35.RData")

# Create ENMTools species objects
K1 <- enmtools.species(species.name = "K1", presence.points = K1_thinned) # presence.points = read.csv("niche-assessment/K1thinned.csv"))
    K1$range <- K1_bias # or `raster(niche-assessment/K1_masked.asc)`
    K1$background.points <- faux_inter_12$points[faux_inter_12$points$species == "dominicensis",]
    check.species(K1)

K2 <- enmtools.species(species.name = "K2", presence.points = K2_thinned) # presence.points = read.csv("niche-assessment/K2thinned.csv"))
    K2$range <- K2_bias
    K2$background.points <- faux_inter_12$points[faux_inter_12$points$species == "South",]
    check.species(K2)

K3 <- enmtools.species(species.name = "K3", presence.points = K3_thinned) # presence.points = read.csv("niche-assessment/K3thinned.csv"))
    K3$range <- K3_bias
    K3$background.points <- faux_inter_13$points[faux_inter_13$points$species == "ignigularis",]
    check.species(K3)

K4 <- enmtools.species(species.name = "K4", presence.points = K4_thinned) # presence.points = read.csv("niche-assessment/K4thinned.csv"))
    K4$range <- K4_bias
    K4$background.points <- faux_inter_24$points[faux_inter_24$points$species == "ravitergum",]
    check.species(K4)

K5 <- enmtools.species(species.name = "K5", presence.points = K5_thinned) # presence.points = read.csv("niche-assessment/K5thinned.csv"))
    K5$range <- K5_bias
    K5$background.points <- faux_inter_35$points[faux_inter_35$points$species == "properus",]
    check.species(K5)

# Perform niche equivalency/identity tests from Warren et al. (2008)
    ## The identity test takes presence points for two species and randomly reassigns them
    ## to each species, builds ENMs for these randomized occurrences, and, through many reps,
    ## estimates a probability distribution for ENM overlap between species under null hypothesis
    ## that species' occurrences are randomly drawn from same distribution.

K_12_id <- identity.test(species.1 = K1, species.2 = K2, nreps = 100, env = env2, type = "mx")
    pdf("niche-assessment/K_1_2.pdf")
        plot(K_12_id)
    dev.off()
    save(K_12_id, file = "niche-assessment/K_1_2_ID.RData")
K_13_id <- identity.test(species.1 = K1, species.2 = K3, nreps = 100, env = env2, type = "mx")
    pdf("niche-assessment/K_1_3.pdf")
        plot(K_13_id)
    dev.off()
    save(K_13_id, file = "niche-assessment/K_1_3_ID.RData")
K_24_id <- identity.test(species.1 = K2, species.2 = K4, nreps = 100, env = env2, type = "mx")
    pdf("niche-assessment/K_2_4.pdf")
        plot(K_24_id)
    dev.off()
    save(K_24_id, file = "niche-assessment/K_2_4_ID.RData")
K_34_id <- identity.test(species.1 = K3, species.2 = K4, nreps = 100, env = env2, type = "mx")
    pdf("niche-assessment/K_3_4.pdf")
        plot(K_34_id)
    dev.off()
    save(K_34_id, file = "niche-assessment/K_3_4_ID.RData")
K_35_id <- identity.test(species.1 = K3, species.2 = K5, nreps = 100, env = env2, type = "mx")
    pdf("niche-assessment/K_3_5.pdf")
        plot(K_35_id)
    dev.off()
    save(K_35_id, file = "niche-assessment/K_3_5_ID.RData")
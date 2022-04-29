
# Load libraries
library(spocc)
library(sp)
library(raster)
library(rgeos)
library(sf)
library(ENMTools)
library(usdm)
library(RStoolbox)
library(tidyverse)

# Load environmental variables
chelsa_clim <- raster::stack(list.files(path = "data/chelsa_new/", pattern = ".asc", full.names = TRUE))
modis_vi <- raster::stack(list.files(path = "data/MODIS_new/", pattern = ".asc", full.names = TRUE))
## Raster stacks have different extents
    ## Since the extent of modis_vi fits within chelsa_clim's extent, crop chelsa_clim
    chelsa_clim <- crop(chelsa_clim, modis_vi@extent)
    chelsa_clim@extent <- raster::alignExtent(chelsa_clim@extent, modis_vi)

## Now, load elevational data that has already been adjusted to the MODIS variable extent
elev_srtm <- raster("data/elevation_new/SRTM_elevation_1km.asc")

## Merge all environmental data layers into single raster stack
env <- raster::stack(elev_srtm, chelsa_clim, modis_vi, full.names = TRUE)
    ## Check rasters for same extent and resolution &
    ## set values to NA if NAs present in any layer
    env <- ENMTools::check.env(env)

env2 <- env[[c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_05",
            "CHELSA_bio10_15", "CHELSA_bio10_16",
            "march_EVI_mean", "may_NDVI_mean")]]

# Create ENMTools species objects
for (k in c(1:3, 5:9)){
    species <- enmtools.species(species.name = paste0("K", k),
                            presence.points = read.csv(paste0("niche-assessment/K", k, "_total_pts.csv")))
    species$range <- background.raster.buffer(species$presence.points, 50000, mask = env2)
    species$background.points <- background.points.buffer(points = species$presence.points,
                                    radius = 1000, n = 1000, mask = env2[[1]])
    print(check.species(species))
        assign(species$species.name, species)
}

# Perform niche equivalency testing
    ## The identity test takes presence points for two species and randomly reassigns them
    ## to each species, builds ENMs for these randomized occurrences, and, through many reps,
    ## estimates a probability distribution for ENM overlap between species under null hypothesis
    ## that species' occurrences are randomly drawn from same distribution.

# A. brevirostris cluster (8) vs A. d. ravitergum cluster (5)
K8_5 <- identity.test(species.1 = K8, species.2 = K5, env = env2, type = "glm", nreps = 500)
    print(K8_7)
    pdf("niche-assessment/K8_5.pdf")
    plot(K8_5)
    dev.off()

# A. brevirostris cluster (8) vs South paleo-island cluster (6)
K8_6 <- identity.test(species.1 = K8, species.2 = K6, env = env2, type = "glm", nreps = 500)
    print(K8_6)
    pdf("niche-assessment/K8_6.pdf")
    plot(K8_6)
    dev.off()

# A. brevirostris cluster (8) vs A. d. ignigularis cluster (7)
K8_7 <- identity.test(species.1 = K8, species.2 = K7, env = env2, type = "glm", nreps = 500)
    print(K8_7)
    pdf("niche-assessment/K8_7.pdf")
    plot(K8_7)
    dev.off()

# A. d. properus cluster (1) vs A. d. ignigularis cluster (7)
# K1_7 <- identity.test(species.1 = K1, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K1_7)
#     pdf("niche-assessment/K1_7.pdf")
#     plot(K1_7)
#     dev.off()

# # A. d. dominicensis II cluster (2) vs A. d. ignigularis cluster (7)
# K2_7 <- identity.test(species.1 = K2, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K2_7)
#     pdf("niche-assessment/K2_7.pdf")
#     plot(K2_7)
#     dev.off()

# # A. d. dominicensis II cluster (2) vs A. d. dominicensis I & IV cluster (3)
# K2_3 <- identity.test(species.1 = K2, species.2 = K3, env = env2, type = "glm", nreps = 500)
#     print(K2_3)
#     pdf("niche-assessment/K2_3.pdf")
#     plot(K2_3)
#     dev.off()

# # A. d. ravitergum cluster (5) vs. A. d. ignigularis cluster (7) 
# K5_7 <- identity.test(species.1 = K5, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K5_7)
#     pdf("niche-assessment/K5_7.pdf")
#     plot(K5_7)
#     dev.off()

# # A. d. ravitergum cluster (5) vs. South paleo-island subsp. cluster (6)
# K5_6 <- identity.test(species.1 = K5, species.2 = K6, env = env2, type = "glm", nreps = 500)
#     print(K5_6)
#     pdf("niche-assessment/K5_6.pdf")
#     plot(K5_6)
#     dev.off()

# # A. d. dominicensis I & IV cluster (3) vs vs. South paleo-island subsp. cluster (6) 
# K3_6 <- identity.test(species.1 = K3, species.2 = K6, env = env2, type = "glm", nreps = 500)
#     print(K3_6)
#     pdf("niche-assessment/K3_6.pdf")
#     plot(K3_6)
#     dev.off()

# # Perform rangebreak tests

# # A. d. properus cluster (1) vs A. d. ignigularis cluster (7)
# K1_7.rbl <- rangebreak.linear(species.1 = K1, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K1_7.rbl)
#     pdf("niche-assessment/K1_7-rangebreak-linear.pdf")
#     plot(K1_7.rbl)
#     dev.off()

# K1_7.rbb <- rangebreak.blob(species.1 = K1, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K1_7.rbb)
#     pdf("niche-assessment/K1_7-rangebreak-blob.pdf")
#     plot(K1_7.rbb)
#     dev.off()

# # A. d. dominicensis II cluster (2) vs A. d. ignigularis cluster (7)
# K2_7.rbl <- rangebreak.linear(species.1 = K2, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K2_7.rbl)
#     pdf("niche-assessment/K2_7-rangebreak-linear.pdf")
#     plot(K2_7.rbl)
#     dev.off()

# K2_7.rbb <- rangebreak.blob(species.1 = K2, species.2 = K7, env = env2, type = "glm", nreps = 500)
#     print(K2_7.rbb)
#     pdf("niche-assessment/K2_7-rangebreak-blob.pdf")
#     plot(K2_7.rbb)
#     dev.off()

# A. d. ravitergum cluster (5) vs. A. d. ignigularis cluster (7)
K5_7.rbl <- rangebreak.linear(species.1 = K5, species.2 = K7, env = env2, type = "glm", nreps = 500)
    print(K5_7.rbl)
    pdf("niche-assessment/K5_7-rangebreak-linear.pdf")
    plot(K5_7.rbl)
    dev.off()

K5_7.rbb <- rangebreak.blob(species.1 = K5, species.2 = K7, env = env2, type = "glm", nreps = 500)
    print(K5_7.rbb)
    pdf("niche-assessment/K5_7-rangebreak-blob.pdf")
    plot(K5_7.rbb)
    dev.off()

# A. d. ravitergum cluster (5) vs. South paleo-island subsp. cluster (6)
K5_6.rbl <- rangebreak.linear(species.1 = K5, species.2 = K6, env = env2, type = "glm", nreps = 500)
    print(K5_6.rbl)
    pdf("niche-assessment/K5_6-rangebreak-linear.pdf")
    plot(K5_6.rbl)
    dev.off()

K5_6.rbb <- rangebreak.blob(species.1 = K5, species.2 = K6, env = env2, type = "glm", nreps = 500)
    print(K5_6.rbb)
    pdf("niche-assessment/K5_6-rangebreak-blob.pdf")
    plot(K5_6.rbb)
    dev.off()

# A. d. dominicensis I & IV cluster (3) vs vs. South paleo-island subsp. cluster (6) 
K3_6.rbl <- rangebreak.linear(species.1 = K3, species.2 = K6, env = env2, type = "glm", nreps = 500)
    print(K3_6.rbl)
    pdf("niche-assessment/K3_6-rangebreak-linear.pdf")
    plot(K3_6.rbl)
    dev.off()

K3_6.rbb <- rangebreak.blob(species.1 = K3, species.2 = K6, env = env2, type = "glm", nreps = 500)
    print(K3_6.rbb)
    pdf("niche-assessment/K3_6-rangebreak-blob.pdf")
    plot(K3_6.rbb)
    dev.off()
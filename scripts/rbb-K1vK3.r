library(tidyverse)
library(raster)
library(rgdal)
library(rgeos)
library(ENMTools)

setwd("/home/tcm0036/distichus-spatial-variation")

# Load environmental variables
env <- raster::stack(list.files(path = "data/env_imputed", pattern = ".asc", full.names = TRUE))    
env2 <- env[[c("CHELSA_bio10_03", "CHELSA_bio10_04",
            "CHELSA_bio10_15", "CHELSA_bio10_16",
            "march_EVI_mean", "may_NDVI_mean")]]
    ENMTools::check.env(env2)

# Create ENMTools species objects
K1 <- enmtools.species(species.name = "K1", presence.points = read.csv("niche-assessment/K_1_thinned.csv"))
    K1$range <- background.raster.buffer(K1$presence.points, radius = 50000, mask = env2[[1]])
    K1$background.points <- background.buffer(points = K1$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K1)

K3 <- enmtools.species(species.name = "K3", presence.points = read.csv("niche-assessment/K_3_thinned.csv"))
    K3$range <- background.raster.buffer(K3$presence.points, radius = 50000, mask = env2[[1]])
    K3$background.points <-background.buffer(points = K3$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K3)

####################################################################################################
# Run rangebreak tests
## Eastern A. d. dominicensis vs A. d. ignigularis
K1K3.rbb <- rangebreak.blob(species.1 = K1, species.2 = K3, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_K1vK3/")
    png("niche-assessment/K_1_3_rbb.png")
        plot(K1K3.rbb)
    dev.off()
    print(summary(K1K3.rbb))
    save(K1K3.rbb, file = "niche-assessment/K1K3.rbb.RData")

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

K2 <- enmtools.species(species.name = "K2", presence.points = read.csv("niche-assessment/K_2_thinned.csv"))
    K2$range <- background.raster.buffer(K2$presence.points, radius = 50000, mask = env2[[1]])
    K2$background.points <- background.buffer(points = K2$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K2)

####################################################################################################
# Run rangebreak tests

# Eastern A. d. dominicensis vs western A. d. dominicensis and South paleo-island
K1K2.rbb <- rangebreak.blob(species.1 = K1, species.2 = K2, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_1_2_rbb.png")
        plot(K1K2.rbb)
    dev.off()
    print(summary(K1K2.rbb))
    save(K1K2.rbb, file = "niche-assessment/K1K2.rbb.RData") 


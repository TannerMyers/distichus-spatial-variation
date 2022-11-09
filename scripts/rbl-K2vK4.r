
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
K2 <- enmtools.species(species.name = "K2", presence.points = read.csv("niche-assessment/K_2_thinned.csv"))
    K2$range <- background.raster.buffer(K2$presence.points, radius = 50000, mask = env2[[1]])
    K2$background.points <- background.buffer(points = K2$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K2)

K4 <- enmtools.species(species.name = "K4", presence.points = read.csv("niche-assessment/K_4_thinned.csv"))
    K4$range <- background.raster.buffer(K4$presence.points, radius = 50000, mask = env2[[1]])
    K4$background.points <- background.buffer(points = K4$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 1000, mask = env2[[1]])
    check.species(K4)

####################################################################################################
# Run rangebreak tests
## A. d. ravitergum vs western A. d. dominicensis and South paleo-island
K4K2.rbl <- rangebreak.linear(species.1 = K4, species.2 = K2, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_K4vK2/")
    png("niche-assessment/K_4_2_rbl.png")
        plot(K4K2.rbl)
    dev.off()
    print(summary(K4K2.rbl))
    save(K4K2.rbl, file = "niche-assessment/K4K2.rbl.RData")

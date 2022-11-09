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

K3 <- enmtools.species(species.name = "K3", presence.points = read.csv("niche-assessment/K_3_thinned.csv"))
    K3$range <- background.raster.buffer(K3$presence.points, radius = 50000, mask = env2[[1]])
    K3$background.points <-background.buffer(points = K3$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K3) 

K5 <- enmtools.species(species.name = "K5", presence.points = read.csv("niche-assessment/K_5_thinned.csv"))
    K5$range <- background.raster.buffer(K5$presence.points, radius = 50000, mask = env2[[1]])
    K5$background.points <- background.buffer(points = K5$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 1000, mask = env2[[1]])
    check.species(K5)

####################################################################################################
# Run rangebreak tests

## A. d. ignigularis vs A. d. properus
K3K5.rbl <- rangebreak.linear(species.1 = K3, species.2 = K5, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_K3vK5/")
    png("niche-assessment/K_3_5_rbl.png")
        plot(K3K5.rbl)
    dev.off()
    print(summary(K3K5.rbl))
    save(K3K5.rbl, file = "niche-assessment/K3K5.rbl.RData")

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
K3K5.rbb <- rangebreak.blob(species.1 = K3, species.2 = K5, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_3_5_rbb.png")
        plot(K3K5.rbb)
    dev.off()
    print(summary(K3K5.rbb))
    save(K3K5.rbb, file = "niche-assessment/K3K5.rbb.RData")

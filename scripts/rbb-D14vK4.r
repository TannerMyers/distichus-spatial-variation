
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
K4 <- enmtools.species(species.name = "K4", presence.points = read.csv("niche-assessment/K_4_thinned.csv"))
    K4$range <- background.raster.buffer(K4$presence.points, radius = 50000, mask = env2[[1]])
    K4$background.points <- background.buffer(points = K4$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 1000, mask = env2[[1]])
    check.species(K4)

d14 <- enmtools.species(species.name = "d14", presence.points = read.csv("niche-assessment/dom14_thinned.csv"))
    d14$range <- background.raster.buffer(d14$presence.points, radius = 50000, mask = env2[[1]])
    d14$background.points <- background.buffer(points = d14$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2)
    check.species(d14)

####################################################################################################
# Run rangebreak tests
K4d14.rbb <- rangebreak.blob(species.1 = K4, species.2 = d14, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_K4vsD14/")
    png("niche-assessment/K_4_d14_rbb.png")
        plot(K4d14.rbb)
    dev.off()
    print(summary(K4d14.rbb))
    save(K4d14.rbb, file = "niche-assessment/K4d14.rbb.RData")
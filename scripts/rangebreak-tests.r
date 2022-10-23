
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

K3 <- enmtools.species(species.name = "K3", presence.points = read.csv("niche-assessment/K_3_thinned.csv"))
    K3$range <- background.raster.buffer(K3$presence.points, radius = 50000, mask = env2[[1]])
    K3$background.points <-background.buffer(points = K3$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(K3)

K4 <- enmtools.species(species.name = "K4", presence.points = read.csv("niche-assessment/K_4_thinned.csv"))
    K4$range <- background.raster.buffer(K4$presence.points, radius = 50000, mask = env2[[1]])
    K4$background.points <- background.buffer(points = K4$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 1000, mask = env2[[1]])
    check.species(K4)

K5 <- enmtools.species(species.name = "K5", presence.points = read.csv("niche-assessment/K_5_thinned.csv"))
    K5$range <- background.raster.buffer(K5$presence.points, radius = 50000, mask = env2[[1]])
    K5$background.points <- background.buffer(points = K5$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 1000, mask = env2[[1]])
    check.species(K5)

# Eastern A. d. dominicensis vs western A. d. dominicensis and South paleo-island
K1K2.rbb <- rangebreak.blob(species.1 = K1, species.2 = K2, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_1_2_rbb.png")
        plot(K1K2.rbb)
    dev.off()
    print(summary(K1K2.rbb))
    save(K1K2.rbb, file = "niche-assessment/K1K2.rbb")

## Eastern A. d. dominicensis vs A. d. ignigularis
K1K3.rbb <- rangebreak.blob(species.1 = K1, species.2 = K3, env = env2, type = "mx", nreps = 1000) 
    png("niche-assessment/K_1_3_rbb.png")
        plot(K1K3.rbb)
    dev.off()
    print(summary(K1K3.rbb))
    save(K1K3.rbb, file = "niche-assessment/K1K3.rbb")

# ## A. d. ignigularis vs A. d. properus
K3K5.rbb <- rangebreak.blob(species.1 = K3, species.2 = K5, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_3_5_rbb.png")
        plot(K3K5.rbb)
    dev.off()
    print(summary(K3K5.rbb))
    save(K3K5.rbb, file = "niche-assessment/K3K5.rbb")

# ## A. d. ignigularis vs A. d. ravitergum
K3K4.rbb <- rangebreak.blob(species.1 = K3, species.2 = K4, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_3_4_rbb.png")
        plot(K3K4.rbb)
    dev.off()
    print(summary(K3K4.rbb))
    save(K3K4.rbb, file = "niche-assessment/K3K4.rbb")

## A. d. ravitergum vs western A. d. dominicensis and South paleo-island
K4K2.rbb <- rangebreak.blob(species.1 = K4, species.2 = K2, env = env2, type = "mx", nreps = 1000)
    png("niche-assessment/K_4_2_rbb.png")
        plot(K4K2.rbb)
    dev.off()
    print(summary(K4K2.rbb))
    save(K4K2.rbb, file = "niche-assessment/K4K2.rbb")

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

spi <- enmtools.species(species.name = "spi", presence.points = read.csv("niche-assessment/south_paleo_thinned.csv"))
    spi$range <- background.raster.buffer(spi$presence.points, radius = 50000, mask = env2[[1]])
    spi$background.points <- background.buffer(points = spi$presence.points[, 1:2], buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(spi)

####################################################################################################
# Run rangebreak tests
K4spi.rbl <- rangebreak.linear(species.1 = K4, species.2 = spi, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_K4vsSPI/linear")
    png("niche-assessment/K_4_spi_rbl.png")
        plot(K4spi.rbl)
    dev.off()
    print(summary(K4spi.rbl))
    save(K4spi.rbl, file = "niche-assessment/K4spi.rbl.RData")


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
spi <- enmtools.species(species.name = "spi", presence.points = read.csv("niche-assessment/south_paleo_thinned.csv"))
    spi$range <- background.raster.buffer(spi$presence.points, radius = 50000, mask = env2[[1]])
    spi$background.points <- background.buffer(points = spi$presence.points[, 1:2], buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2[[1]])
    check.species(spi)

d14 <- enmtools.species(species.name = "d14", presence.points = read.csv("niche-assessment/dom14_thinned.csv"))
    d14$range <- background.raster.buffer(d14$presence.points, radius = 50000, mask = env2[[1]])
    d14$background.points <- background.buffer(points = d14$presence.points, buffer.width = 20000, buffer.type = "convhull", return.type = "points", n = 3000, mask = env2)
    check.species(d14)

####################################################################################################
# Run rangebreak tests
spid14.rbl <- rangebreak.linear(species.1 = spi, species.2 = d14, env = env2, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/niche/rangebreak_SPIvsD14/")
    png("niche-assessment/K_4_spi_rbl.png")
        plot(spid14.rbl)
    dev.off()
    print(summary(spid14.rbl))
    save(spid14.rbl, file = "niche-assessment/spid14.rbl.RData")

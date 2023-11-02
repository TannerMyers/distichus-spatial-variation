library(tidyverse)
library(terra)
library(sf)
library(sp)
library(ENMTools)

setwd("/home/tcm0036/distichus-spatial-variation")

# Load environmental variables
env <- terra::rast(list.files(path = "data/selected-vars/imputed/", pattern = ".asc$", full.names = TRUE))
# crs(env) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
names(env) <- c("bio18", "bio3", "bio4", "bio9", "clt_min", "cmi_range", "elevation")
    ENMTools::check.env(env)

## eastern dominicensis
K2 <- enmtools.species(species.name = "K2", presence.points = vect(read.csv("niche-assessment/K2_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K2$presence.points) <- crs(env)
    K2$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/K2_alpha_poly.shp"))
    K2$background.points <- background.points.buffer(points = K2$presence.points, radius = 15000, n = 3000, mask = K2$range)
    check.species(K2) 

## ignigularis
K5 <- enmtools.species(species.name = "K5", presence.points = vect(read.csv("niche-assessment/K5_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K5$presence.points) <- crs(env)
    K5$range <- terra::mask(x = env[[1]], mask = vect("niche-assessment/ranges/K5_alpha_poly.shp"))
    K5$background.points <- background.points.buffer(points = K5$presence.points, radius = 15000, n = 3000, mask = K5$range)
    check.species(K5)

####################################################################################################
# Run rangebreak tests
## A. d. ignigularis vs East A.d. dominicensis
dom_ign <- rangebreak.blob(species.1 = K2, species.2 = K5, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-east-dom-ign/")
    png("niche-assessment/dom_ign_rbb.png")
        plot(dom_ign)
    dev.off()
    print(summary(dom_ign))
    save(dom_ign, file = "niche-assessment/dom_ign_blob.RData")

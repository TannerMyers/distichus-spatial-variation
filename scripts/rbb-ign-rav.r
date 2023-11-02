
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

## ravitergum 
K3 <- enmtools.species(species.name = "K3", presence.points = vect(read.csv("niche-assessment/K3_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K3$presence.points) <- crs(env)
    K3$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/K3_alpha_poly.shp"))
    K3$background.points <- background.points.buffer(points = K3$presence.points, radius = 15000, n = 3000, mask = K3$range)
    check.species(K3)

## ignigularis 
K5 <- enmtools.species(species.name = "K5", presence.points = vect(read.csv("niche-assessment/K5_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K5$presence.points) <- crs(env)
    K5$range <- terra::mask(x = env[[1]], mask = vect("niche-assessment/ranges/K5_alpha_poly.shp"))
    K5$background.points <- background.points.buffer(points = K5$presence.points, radius = 15000, n = 3000, mask = K5$range)
    check.species(K5)

####################################################################################################
## Run range-break test
## A. d. ignigularis vs A. d. ravitergum
ign_rav.rbb <- rangebreak.blob(species.1 = K3, species.2 = K5, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-ign-rav/")
    png("niche-assessment/ign_rav_rbb.png")
        plot(ign_rav.rbb)
    dev.off()
    print(summary(ign_rav.rbb))
    save(ign_rav.rbb, file = "niche-assessment/ign_rav_blob.RData")


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

# Create ENMTools species objects

## western dominicensis
west_dom <- enmtools.species(species.name = "western_dominicensis", presence.points = vect(read.csv("niche-assessment/dom14_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(west_dom$presence.points) <- crs(env)
    west_dom$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/dom14_poly.shp"))
    west_dom$background.points <- background.points.buffer(points = west_dom$presence.points, radius = 22000, n = 3000, mask = west_dom$range)
    check.species(west_dom) 

## eastern dominicensis
K2 <- enmtools.species(species.name = "K2", presence.points = vect(read.csv("niche-assessment/K2_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K2$presence.points) <- crs(env)
    K2$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/K2_alpha_poly.shp"))
    K2$background.points <- background.points.buffer(points = K2$presence.points, radius = 15000, n = 3000, mask = K2$range)
    check.species(K2)

####################################################################################################
## Run range-break test
## East A. d. dominicensis vs West A. d. dominicensis
ew_dom <- rangebreak.blob(species.1 = west_dom, species.2 = K2, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-east-west-dom/")
    png("niche-assessment/ew-dom-rbb.png")
        plot(ew_dom)
    dev.off()
    print(summary(ew_dom))
    save(ew_dom, file = "niche-assessment/ew_dom_blob.RData")

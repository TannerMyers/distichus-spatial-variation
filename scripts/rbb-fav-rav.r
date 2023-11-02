
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

## ravitergum 
K3 <- enmtools.species(species.name = "K3", presence.points = vect(read.csv("niche-assessment/K3_alpha_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(K3$presence.points) <- crs(env)
    K3$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/K3_alpha_poly.shp"))
    K3$background.points <- background.points.buffer(points = K3$presence.points, radius = 15000, n = 3000, mask = K3$range)
    check.species(K3)

fav <- enmtools.species(species.name = "favillarum", presence.points = vect(read.csv("niche-assessment/fav_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(fav$presence.points) <- crs(env)
    fav$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/fav_poly.shp"))
    fav$background.points <- background.points.buffer(points = fav$presence.points, radius = 2000, n = 3000, mask = fav$range)
    check.species(fav)

####################################################################################################
## Run range-break test
## South paleo-island vs A. d. ravitergum
rav_fav <- rangebreak.blob(species.1 = K3, species.2 = fav, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-rav-fav/")
    png("niche-assessment/rav-fav-rbb.png")
        plot(rav_fav)
    dev.off()
    print(summary(rav_fav))
    save(rav_fav, file = "niche-assessment/rav_fav_blob.RData")

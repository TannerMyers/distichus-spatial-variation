
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
west_dom <- enmtools.species(species.name = "western_dominicensis", presence.points = vect(read.csv("niche-assessment/dom14_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(west_dom$presence.points) <- crs(env)
    west_dom$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/dom14_poly.shp"))
    west_dom$background.points <- background.points.buffer(points = west_dom$presence.points, radius = 22000, n = 3000, mask = west_dom$range)
    check.species(west_dom)   

tiburon <- enmtools.species(species.name = "tiburon", presence.points = vect(read.csv("niche-assessment/tiburon_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(tiburon$presence.points) <- crs(env)
    tiburon$range <- background.raster.buffer(tiburon$presence.points, radius = 33000, mask = env[[1]]) 
    tiburon$background.points <- background.points.buffer(points = tiburon$presence.points, radius = 32000, n = 3000, mask = tiburon$range)
    check.species(tiburon)

####################################################################################################
## Run range-break test
## A. d. favillarum vs West A. d. dominicensis
w_dom_tiburon <- rangebreak.blob(species.1 = west_dom, species.2 = tiburon, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-west-dom-tiburon/")
    png("niche-assessment/w-dom-sur-rbb.png")
        plot(w_dom_tiburon)
    dev.off()
    print(summary(w_dom_tiburon))
    save(w_dom_tiburon, file = "niche-assessment/w_dom_tiburon_blob.RData")

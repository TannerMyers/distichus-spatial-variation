
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

fav <- enmtools.species(species.name = "favillarum", presence.points = vect(read.csv("niche-assessment/fav_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(fav$presence.points) <- crs(env)
    fav$range <- terra::mask(x = env[[1]],  mask = vect("niche-assessment/ranges/fav_poly.shp"))
    fav$background.points <- background.points.buffer(points = fav$presence.points, radius = 2000, n = 3000, mask = fav$range)
    check.species(fav)

tiburon <- enmtools.species(species.name = "tiburon", presence.points = vect(read.csv("niche-assessment/tiburon_poly-thinned-points.csv"), geom = c("longitude", "latitude")))
    crs(tiburon$presence.points) <- crs(env)
    tiburon$range <- background.raster.buffer(tiburon$presence.points, radius = 33000, mask = env[[1]]) 
    tiburon$background.points <- background.points.buffer(points = tiburon$presence.points, radius = 32000, n = 3000, mask = tiburon$range)
    check.species(tiburon)

####################################################################################################
## Run range-break test
## Tiburón peninsula endemic subspecies vs A. d. favillarum
fav_tiburon <- rangebreak.blob(species.1 = fav, species.2 = tiburon, env = env, type = "mx", nreps = 1000,
                            bg.source = "points", low.memory = TRUE, verbose = TRUE,
                            rep.dir = "/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/niche/rangebreak-fav-tiburon/")
    png("niche-assessment/fav-tiburon-rbb.png")
        plot(fav_tiburon)
    dev.off()
    print(summary(fav_tiburon))
    save(fav_tiburon, file = "niche-assessment/fav_tiburon_blob.RData")

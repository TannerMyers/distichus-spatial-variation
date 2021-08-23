#
# Process MERRAclim data for Hispaniola
# From: https://rdrr.io/github/RS-eco/bdc/src/data-raw/merraclim.R

setwd("~/Dropbox/Distichus_Project/distichus_spatial_morphological_variation/")

# Load libraries
library(raster)
library(dplyr)

# Set file directory
  filedir <- "~/Dropbox/Distichus_Project/distichus_spatial_morphological_variation/environmental-data/MERRAclim/"

# Create file names
   df <- expand.grid(vars=c("min", "max", "mean"), years=c("00s","90s","80s"),
                 res=c("2_5m"#, "5m", "10m"
                        ))
   df$files <- paste0(df$res, "_", df$vars, "_", df$years, ".zip")

# Unzip files
   #lapply(df$files, function(file){unzip(paste0(filedir, file), exdir = filedir)})

# Load files
  dat <- lapply(df$files, function(file){raster::stack(
              list.files(filedir, pattern=paste0("^", sub(".zip", "", file)), 
                         full.names=T)[1:15])})

# Crop by extent of Hispaniola
coords <- read.csv('coords.csv', header=TRUE) # Loading dataframe for all 110 localities with measured distichus
min.lat <- floor(min(coords$Latitude))
max.lat <- ceiling(max(coords$Latitude))
max.lon <- ceiling(max(coords$Longitude))
min.lon <- floor(min(coords$Longitude))
geog.extent <- raster::extent(x = c(min.lon, max.lon, min.lat, max.lat))

dat_Dom <- lapply(dat, function(x) raster::crop(x, geog.extent, snap="out"))

# Skipping a bunch of lines because I don't understand them and don't have a grid
# file to resample rasters with

merraclim_Dom <- lapply(1:length(dat_Dom), function(x){
  dat <- as.data.frame(raster::rasterToPoints(dat_Dom[[x]]))
  dat$var <- df$vars[[x]]
  dat$year <- df$years[[x]]
  dat$res <- df$res[[x]]
  colnames(dat) <- sub(paste0("X", df$res[[x]], "_", df$vars[[x]], "_", df$years[[x]], "_"), "", colnames(dat))
  return(dat)
})

merraclim_Dom <- dplyr::bind_rows(merraclim_Dom)
colnames(merraclim_Dom)
head(merraclim_Dom)

merraclim_2.5m_Dom <- merraclim_Dom %>% filter(res == "2_5m") %>% select(-res) # repeat if you have 5m or 10m

# Save resulting dataframe of values
save(merraclim_2.5m_Dom, file="environmental-data/MERRAclim/merraclim_2.5m_hispaniola.rda", 
     compress = "xz")

#load("environmental-data/MERRAclim/merraclim_2.5m_hispaniola.rda")
clim <- raster::rasterFromXYZ(merraclim_2.5m_Dom)
  # how to account for year? 
  # Average based on locality? Use x (longitude)
climValues <- raster::extract(clim, coords[,1:2]) 
str(climValues)




###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# Code obtained and modified from https://github.com/townpeterson/vespa/tree/master/Rcode

# This script downloads occurrence records for all of the distichus lineages identified by MacGuigan et al. (2017)
# found on Hispaniola: favillarum, ignigularis, ravitergum, dominicensis2, properus, and the Tiburon penninsula
# endemic subspecies.

###############

# load packages

#library(devtools)
#devtools::install_github("marlonecobos/ellipsenm")
library(spocc)
library(rgbif)
library(rvertnet)
library(maps)
library(ellipsenm)
library(rgdal)
library(raster)
library(sp)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggmap)
library(httpgd)
library(tidyverse)

# Directories
working_dir <- getwd()
setwd(working_dir)

output_dir <- "~/Dropbox/Distichus_Project/distichus-spatial-project/data/GBIF/" # Add where you want occurrences to go
dir.create(output_dir)

#-------------------------------------------------------------------------------

#----------------------------Downloading occurrences----------------------------

# getting occurrences from GBIF
ad <- occ("Anolis distichus", limit = 100000)
gbif_occs <- ad$gbif$data$Anolis_distichus
head(gbif_occs)

## write gbif_occs to file to create a derived dataset – this gives us a citable DOI for the GBIF data we're using
write_csv(gbif_occs, "data/GBIF/Anolis_distichus_GBIF_download.csv")

sp_search <- occ_search(taxonKey = ad$gbif$Anolis_distichus$taxonKey[1])
cit <- gbif_citation(sp_search)
sink('data/GBIF/gbif_ref.txt')
sapply(cit, print)
sink()

## saving initial data with various columns from GBIF
colnames(gbif_occs)
columns <- c("name", "longitude", "latitude", "issues", "scientificName", 
             "coordinateUncertaintyInMeters", "year", "month", "day", 
             "countryCode", "locality", "elevation")

occs <- gbif_occs[, columns]

# getting occurrences from VertNet

# Plot occurrence points
hispaniola <- map_data("world")%>%subset(region == c("Haiti", "Dominican Republic"))
ggplot() +
    geom_polygon(data=hispaniola, aes(x=long, y=lat, group=group), color='lightgrey', fill=NA) +
    geom_point(data=occs, aes(x=longitude, y=latitude), shape=1, size=2, alpha=0.85, color="green") + 
    coord_fixed(xlim=c(-75,-67.5), ylim=c(17,20), ratio=1.3) + theme_nothing()

# Write csv containing occurrence data for Anolis distichus 
write.csv(occs, "data/GBIF/occurrences_distichus_GBIF.csv", row.names = FALSE) # Maybe split script here. There's a lot of code after this and there's no need to repeat or alter the steps before

#---------------------------------------------------------------------------------

#--------------------- Data cleaning ---------------------------------------------
occs <- read.csv("data/GBIF/occurrences_distichus_GBIF.csv", header = TRUE)
colnames(occs)

## subset columns of interest
occ <- occs[, c("scientificName", "longitude", "latitude", "year")]
occ <- na.omit(occ) # exclude NAs here after reducing number of fields in dataframe

## excluding 0, 0 coordinates 
occ <- occ[occ$longitude != 0 & occ$latitude != 0, ]

## excluding duplicates
occ <- occ[!duplicated(paste(occ$longitude, occ$latitude)), ]

hgd()
hgd_browse()
maps::map()
points(occ[, 2:3], col = "red", pch = 19)
axis(side = 2)
axis(side =1)

# excluding records outside Caribbean
occ <- occ[occ$longitude < 0,]
occ <- occ[occ$latitude < 30,]

# Omit samples not from Hispaniola
occ <- occ[occ$latitude>16.5,]
occ <- occ[occ$latitude<21,]
occ <- occ[occ$longitude<(-68),]
occ<- occ[occ$longitude>(-80),]

maps::map(xlim=c(-100,-60),ylim=c(-0,30), interior=TRUE)
points(occ[, 2:3], col = "red", pch = 19)
axis(side=2)
axis(side=1)

write.csv(occ, paste0(output_dir, "A_distichus_clean.csv"), row.names=FALSE)

#-------------------------------------------------------------------------------

#-----------------------------Thinning------------------------------------------

## make directory to hold input files for ecological niche models
## In this case, thinned occurrence records
dir.create("enm/thinned-datasets")

## spatially thin occurrences by 30 km
occ_thin <- thin_data(occ, "longitude", "latitude", thin_distance = 30,
                    save =TRUE, name = "enm/thinned-datasets/A_distichus_30km")

maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_thin[,2:3], col = "red", pch =19)

#-------------------------------------------------------------------------------

#-----------------Splitting training and testing data---------------------------
# model calibration folder
dir.create("enm/Model_calibration")
dir.create("enm/Model_calibration/All_Records_30km_thin")

# split distance based thinned data for 5 exercises of model calibration
n <- 1:5
splits <- lapply(n, function(x) {
  set.seed(x)
  occsp_50 <- split_data(occ_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/All_Records_30km_thin/Adist", x))
})

#-------------------------------------------------------------------------------

#-----------------Subset data by genetically distinct lineages -----------------

# I want to curate occurrences from A. d. favillarum, A. d. ignigularis, A. d. ravitergum,
# A. d. properus and A. d. sejunctus, A. d. dominicensis (1 & 2), and the Haitian endemic 
# subspecies.

## Do favillarum subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName==c("Anolis distichus favillarum Schwartz, 1968"),2:3], col = "red", pch =19)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_fav <- occ[occ$latitude<=18.21,]
occ_fav <- occ_fav[occ_fav$latitude>=17.96,]
occ_fav <- occ_fav[occ_fav$longitude<=(-71.07444),]
occ_fav <- occ_fav[occ_fav$longitude>=(-71.29306),] # done

## Do ignigularis subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName=="Anolis distichus ignigularis Mertens, 1939",2:3], col = "red", pch =19)
points(occ[occ$scientificName=="Anolis ignigularis Glor & Laport, 2012",2:3], col = "red", pch =19)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_ig <- occ[occ$latitude<=19.16,]
occ_ig <- occ_ig[occ_ig$latitude>=18.03407,]
occ_ig <- occ_ig[occ_ig$longitude<=(-69),]
occ_ig <- occ_ig[occ_ig$longitude>=(-70.43866),] # done

# Remove occurrence in ocean
occ_ig <- occ_ig %>% filter(latitude != 18.23861)

# For the Samaná penninsula population of A. d. ignigularis
occ_ig_Samana <- occ[occ$longitude>=(-69.66),]
occ_ig_Samana <- occ_ig_Samana[occ_ig_Samana$latitude>=19.16,]

occ_ig <- bind_rows(occ_ig_Samana, occ_ig)

## Do ravitergum subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName==c("Anolis ravitergum Glor & Laport, 2012"), 2:3], 
      col = "red", pch =19)
points(occ[occ$scientificName==c("Anolis distichus ravitergum Schwartz, 1968"), 2:3], 
      col = "red", pch =19)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_rav <- occ[occ$latitude<=18.4692,]
occ_rav <- occ_rav[occ_rav$latitude>=18.2383,]
occ_rav <- occ_rav[occ_rav$longitude<=(-70.32293),]
occ_rav <- occ_rav[occ_rav$longitude>=(-71.58279),] # done

## Do properus/sejunctus subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName=="Anolis distichus sejunctus Schwartz, 1968", 2:3], col = "red", pch =19)
points(occ[occ$scientificName=="Anolis distichus properus Schwartz, 1968", 2:3], col = "red", pch =19) 
points(occ[occ$scientificName=="Anolis properus Glor & Laport, 2012",2:3], col = "red", pch =19)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_prop <- occ[occ$latitude<=18.9,]
occ_prop <- occ[occ$latitude>=18.3,]
occ_prop <- occ_prop[occ_prop$longitude>=(-69),]
occ_prop <- occ_prop[occ_prop$longitude<(-68.2251),]

## Do Haitian subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName=="Anolis distichus vinosus Schwartz 1968",2:3], col = "red", pch =19)
points(occ[occ$scientificName=="Anolis distichus aurifer Schwartz, 1968",2:3], col = "red", pch =19) 
points(occ[occ$scientificName=="Anolis distichus suppar Schwartz, 1968",2:3],col = "red", pch =19)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_haiti <- occ[occ$latitude<18.72,]
occ_haiti <- occ_haiti[occ_haiti$longitude<(-71.8),] # done

## Do dominicensis subset
maps::map(xlim = range(occ$longitude), ylim = range(occ$latitude), interior =TRUE)
points(occ[occ$scientificName==c("Anolis distichus dominicensis Reinhardt And Lütken, 1863"),2:3], col = "red", pch =19) 
points(occ[occ$scientificName==c("Anolis dominicensis Reinhardt And Lütken, 1863"),2:3], col = "blue", pch =20)
points(occ[occ$scientificName==c("Anolis distichus albidogularis Mertens, 1939"),2:3], col = "yellow", pch =21)
points(occ[,2:3], col = "blue", pch =1) # check if other observations are in focal area without label

occ_dom <- occ[occ$latitude>=18.72,]
occ_dom <- occ_dom[occ_dom$latitude<=20,]
occ_dom <- occ_dom[occ_dom$longitude<=(-69.62),]

occ_dom <- dplyr::anti_join(occ_dom, occ_ig, by=c("latitude","longitude")) # eliminate overlapping ignigularis samples

write.csv(occ_fav, paste0(output_dir,"A_d_favillarum.csv"), row.names=FALSE)
write.csv(occ_dom, paste0(output_dir,"A_d_dominicensis_clean.csv"), row.names=FALSE)
write.csv(occ_haiti, paste0(output_dir,"A_distichus_Tiburon.csv"), row.names=FALSE)
write.csv(occ_prop, paste0(output_dir,"A_d_properus.csv"), row.names=FALSE)
write.csv(occ_rav, paste0(output_dir,"A_d_ravitergum.csv"), row.names=FALSE)
write.csv(occ_ig, paste0(output_dir,"A_d_ignigularis.csv"), row.names=FALSE)

#-------------------------------------------------------------------------------

#-----------------------------Thinning------------------------------------------

## spatially thin A. d. ignigularis occurrences by 10 km
occ_ig_thin <- thin_data(occ_ig, "longitude", "latitude", thin_distance = 10,
                    save =TRUE, name = "enm/thinned-datasets/A_d_ignigularis_10km")

maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_ig_thin[,2:3], col = "orange", pch =19)

## spatially thin A. d. ravitergum occurrences by 10 km
occ_rav_thin <- thin_data(occ_rav, "longitude", "latitude", thin_distance = 10,
                    save =TRUE, name = "enm/thinned-datasets/A_d_ravitergum_10km")

#maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_rav_thin[,2:3], col = "yellow", pch =19)

## spatially thin A. d. properus occurrences by 10 km
occ_prop_thin <- thin_data(occ_prop, "longitude", "latitude", thin_distance = 10,
                    save =TRUE, name = "enm/thinned-datasets/A_d_properus_10km")

#maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_prop_thin[,2:3], col = "gray", pch =19)

## spatially thin A. d. favillarum occurrences by 10 km
occ_fav_thin <- thin_data(occ_fav, "longitude", "latitude", thin_distance = 3,
                    save =TRUE, name = "enm/thinned-datasets/A_d_favillarum_3km")

#maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_fav_thin[,2:3], col = "purple", pch =19)


## spatially thin A. d. dominicensis occurrences by 10 km
occ_dom_thin <- thin_data(occ_dom, "longitude", "latitude", thin_distance = 15,
                    save =TRUE, name = "enm/thinned-datasets/A_d_dominicensis_30km")

#maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_dom_thin[,2:3], col = "maroon", pch =19)

## spatially thin Haitian A. distichus subspecies occurrences by 10 km
occ_haiti_thin <- thin_data(occ_haiti, "longitude", "latitude", thin_distance = 10,
                    save =TRUE, name = "enm/thinned-datasets/A_distichus_Tiburon_10km")

#maps::map(xlim = range(occ_thin$longitude), ylim = range(occ_thin$latitude), interior =TRUE)
points(occ_haiti_thin[,2:3], col = "blue", pch =19)

#-------------------------------------------------------------------------------

#-----------------Splitting training and testing data---------------------------
# model calibration folder
dir.create("enm/Model_calibration/fav_records_thin")
dir.create("enm/Model_calibration/rav_records_thin")
dir.create("enm/Model_calibration/ig_records_thin")
dir.create("enm/Model_calibration/dom_records_thin")
dir.create("enm/Model_calibration/prop_records_thin")
dir.create("enm/Model_calibration/haiti_records_thin")

# split distance based thinned data for 5 exercises of model calibration
n <- 1:5
fav_splits <- lapply(n, function(x) {
  set.seed(x)
  occfav_3 <- split_data(occ_fav_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/fav_records_thin/fav", x))
})

n <- 1:5
rav_splits <- lapply(n, function(x) {
  set.seed(x)
  occrav_10 <- split_data(occ_rav_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/rav_records_thin/rav", x))
})

n <- 1:5
ig_splits <- lapply(n, function(x) {
  set.seed(x)
  occig_10 <- split_data(occ_ig_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/ig_records_thin/ig", x))
})

n <- 1:5
prop_splits <- lapply(n, function(x) {
  set.seed(x)
  occprop_10 <- split_data(occ_prop_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/prop_records_thin/prop", x))
})

n <- 1:5
dom_splits <- lapply(n, function(x) {
  set.seed(x)
  occdom_15 <- split_data(occ_dom_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/dom_records_thin/dom", x))
})

n <- 1:5
haiti_splits <- lapply(n, function(x) {
  set.seed(x)
  occhaiti_10 <- split_data(occ_haiti_thin, method = "random", longitude = "longitude", 
                         latitude = "latitude", train_proportion = 0.51, save = T, # 0.51 trick for getting training with 1 record more than testing
                         name = paste0("enm/Model_calibration/haiti_records_thin/haiti", x))
})

#-------------------------------------------------------------------------------
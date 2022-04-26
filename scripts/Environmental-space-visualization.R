###############

# Project: Do environmental and geographic variation explain
#          morphological variation in the Hispaniolan bark anole,
#          Anolis distichus?

# Authors:
# Tanner Myers, Pietro de Mello, and Rich Glor

# This script plots values for the first and twelfth bioclimatic variables (i.e., annual mean temperature
# and precipitation) from occurrence records of the different lineages of Anolis distichus identified by MacGuigan et al. (2017) 
# obtained from GBIF, as well as samples with morphological data analyzed by Myers et al. (2020). This is an
# exploratory analysis to begin testing if the different lineages possess distinct niches and to compare different
# sources of bioclimatic variables (WorldClim, MERRAClim, and CHELSA).

# Load packages
library(tidyverse)

# Directories
working_dir <- getwd()
setwd(working_dir)

# Load data for individuals that have been measured 
data <- read_csv("complete-data.csv")
str(data)
head(data)

     ## Designate population colors for plotting
     data$Color[data$Subspecies=="aurifer"] <-"peru" 
     data$Color[data$Subspecies=="vinosus"] <-"firebrick3" 
     data$Color[data$Subspecies=="suppar"] <-"purple1" 
     data$Color[data$Subspecies=="favillarum"] <-"dodgerblue" 
     data$Color[data$Subspecies=="ravitergum"] <-"yellow" 
     data$Color[data$Subspecies=="dominicensis"] <-"bisque" 
     data$Color[data$Subspecies=="ignigularis"] <-"darkorange2" 
     data$Color[data$Subspecies=="properus"] <-"seagreen3"   

columns <- c("Subspecies", "Longitude", "Latitude", "Color")

df <- data[,columns]

# Load processed MERRAClim bioclimatic variables at 2.5 arc seconds resolution
mc_bioclim <- stack(list.files(path="data/merra_new", pattern = ".asc", full.names = T))

wc_bioclim <- stack(list.files(path="data/wc0.5", pattern = ".bil", full.names = T)) 

chelsa_bioclim <- stack(list.files(path="data/chelsa_new/", pattern=".asc", full.names=TRUE))

# Get values for all raster points to plot as background
    ## First, bio1
     mc_bio1 <- mc_bioclim$bio1
     mc_bio1_points <- rasterToPoints(mc_bio1)
     wc_bio1 <- wc_bioclim$bio1_23
     wc_bio1_points <- rasterToPoints(wc_bio1)
     chelsa_bio1 <- chelsa_bioclim$CHELSA_bio10_01
     chelsa_bio1_points <- rasterToPoints(chelsa_bio1)
    
    ## Now, bio12 
     mc_bio12 <-mc_bioclim$bio12
     mc_bio12_points <- rasterToPoints(mc_bio12)
     wc_bio12 <- wc_bioclim$bio12_23
     wc_bio12_points <- rasterToPoints(wc_bio12)
     chelsa_bio12 <- chelsa_bioclim$CHELSA_bio10_12
     chelsa_bio12_points <- rasterToPoints(chelsa_bio12)

# Extract raster values for the localities in `df`
mc_bc_values <- raster::extract(mc_bioclim,data.frame(df$Longitude, df$Latitude))
wc_bc_values <- raster::extract(wc_bioclim,data.frame[df$Longitude, df$Latitude])
chelsa_bc_values <- raster::extract(chelsa_bioclim, data.frame(df$Longitude, df$Latitude))

# Combine into one dataframe
df <- cbind(df, mc_bc_values, wc_bc_values, chelsa_bc_values)

leg.txt <- c(expression(italic("ignigularis")),
             expression(italic("ravitergum")),
             expression(italic("favillarum")),
             expression(italic("properus")),
             expression(italic("dominicensis")),
             expression(italic("aurifer")),
             expression(italic("suppar")),
             expression(italic("vinosus")))

pdf("CHELSA_environmental_biplot.pdf")
plot(chelsa_bio12_points[,3] ~ chelsa_bio1_points[,3], xlab="BIO1", ylab="BIO12", main="Environmental Space (CHELSA)")
     points(df$CHELSA_bio10_12 ~ df$CHELSA_bio10_01, col=df$Color, pch=16, cex=1.25)
     legend("topleft", leg.txt, col=unique(df$Color), pch=16, cex=1)
     dev.off()

pdf("MERRAclim-environmental-biplot.pdf")
#hgd()
#hgd_browse()
plot(mc_bio12_points[,3] ~ mc_bio1_points[,3], 
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (MERRAClim -- 10m 2000s Mean)")
points(df$bio12 ~ df$bio1, 
      col=df$Color, pch=16, cex=1.25)
legend("topleft", leg.txt, col=unique(df$Color), pch =16,
       cex=1)
dev.off()
  
pdf("Worldclim-environmental-biplot.pdf")
plot(wc_bio12_points[,3] ~ wc_bio1_points[,3],
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (WorldClim)")
points(df$bio12_23 ~ df$bio1_23, 
       col=df$Color, pch=16, cex=1.25)
legend("topleft", leg.txt, col=unique(df$Color), pch =16,
         cex=1)
dev.off()


### Now do GBIF data
dom <- read_csv("data/GBIF/A_d_dominicensis_clean.csv")
ign <- read_csv("data/GBIF/A_d_ignigularis.csv")
fav <- read_csv("data/GBIF/A_d_favillarum.csv")
prop <- read_csv("data/GBIF/A_d_properus.csv")
rav <- read_csv("data/GBIF/A_d_ravitergum.csv")
haiti <- read_csv("data/GBIF/A_distichus_Tiburon.csv")

dom$MacGuigan <- "dominicensis"
ign$MacGuigan <- "ignigularis"
fav$MacGuigan <- "favillarum"
rav$MacGuigan <- "ravitergum"
haiti$MacGuigan <- "Tiburon"
prop$MacGuigan <- "properus"

# Merge into single dataframe
data2 <-rbind(ign, rav, fav, prop, dom, haiti)

## Designate population colors for plotting
     data2$Color[data2$Subspecies=="aurifer"] <-"peru" 
     data2$Color[data2$Subspecies=="vinosus"] <-"firebrick3" 
     data2$Color[data2$Subspecies=="suppar"] <-"purple1" 
     data2$Color[data2$Subspecies=="favillarum"] <-"dodgerblue" 
     data2$Color[data2$Subspecies=="ravitergum"] <-"yellow" 
     data2$Color[data2$Subspecies=="dominicensis"] <-"bisque" 
     data2$Color[data2$Subspecies=="ignigularis"] <-"darkorange2" 
     data2$Color[data2$Subspecies=="properus"] <-"seagreen3"   

mc_bc_values2 <- raster::extract(mc_bioclim,data.frame(data2$longitude, data2$latitude))
wc_bc_values2 <- raster::extract(wc_bioclim,data.frame[data2$longitude, data2$latitude])
chelsa_bc_values2 <- raster::extract(chelsa_bioclim,data.frame(data2$longitude, data2$latitude))

data2 <- cbind(data2, mc_bc_values2, wc_bc_values2, chelsa_bc_values2)

pdf("GBIF-MERRAclim-environmental-biplot.pdf")
#hgd()
#hgd_browse()
plot(mc_bio12_points[,3] ~ mc_bio1_points[,3], 
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (MERRAClim -- 10m 2000s Mean)")
points(data2$bio12 ~ data2$bio1, 
      col=data2$Color, pch=16, cex=1.25)
legend("topleft", leg.txt, col=unique(data2$Color), pch =16,
       cex=1)
dev.off()
  
pdf("GBIF-Worldclim-environmental-biplot.pdf")
plot(wc_bio12_points[,3] ~ wc_bio1_points[,3],
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (WorldClim)")
points(data2$bio12_23 ~ data2$bio1_23, 
       col=df$Color, pch=16, cex=1.25)
legend("topleft", leg.txt, col=unique(data2$Color), pch =16,
         cex=1)
dev.off()

pdf("GBIF-CHELSA-environmental-biplot.pdf")
plot(chelsa_bio12_points[,3] ~ chelsa_bio1_points[,3],
     xlab="Annual Mean Temperature", ylab="Annual Precipitation",
     main="Environmental Space (WorldClim)")
points(data2$CHELSA_bio10_12 ~ data2$CHELSA_bio10_01, 
       col=data2$Color, pch=16, cex=1.25)     
legend("topleft", leg.txt, col=unique(data2$Color), pch =16,
         cex=1)     
dev.off()
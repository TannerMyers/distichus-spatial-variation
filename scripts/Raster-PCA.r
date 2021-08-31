
# directories 
working_dir <- getwd()
setwd(working_dir)

# libraries
library(raster)
library(rgeos)
library(rgdal)
library(RStoolbox)
library(tidyverse)

## load environmental variables as a raster stack 
env <- raster::stack(list.files())
plot(env)

## perform principal component analysis on rasters
pca <- rasterPCA(env, spca=TRUE)
    ## Visualize PC rasters
    plot(pca$map)

## Check loadings & eigenvalues for top 5 loadings
knitr::kable(round(pca$model$loadings[,1:5],5)) # loadings
summary(pca$model) # eigenvalues


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
env <- raster::stack(list.files("enm/calibration-areas/A_distichus-calibration-area", pattern=".asc", full.names=TRUE))
plot(env)

## perform principal component analysis on rasters
pca <- rasterPCA(env, spca=TRUE)
    ## Visualize PC rasters
    plot(pca$map)

## Check loadings & eigenvalues for top 5 loadings
summary(pca$model) # eigenvalues
knitr::kable(round(pca$model$loadings[,1:5],5)) # loadings

dir.create("pcas")

spca1 <- spca(layers_stack=env, layers_to_proj = NULL, # change to rasters of world -- unsure if this is necessary
            sv_dir = "pcas", layers_format = ".asc",
            sv_proj_dir = NULL)

spcs <- readRDS("pcas/pca_object21_10_06_16_17.rds")

spcs2 <- summary(spcs)

png(filename="pcas/screeplot_chelsa.png", width=1200*1.3, height=1200*1.3, res=300)
plot(spcs2$importance[3, 1:5]*100, xlab="Principal Component", ylab="Percentage of variance explained",
    type="b", frame.plot=TRUE, cex=1.5)
points(spcs2$importance[2,1:5]*100, pch=17, cex=1.5)
lines(spcs2$importance[2,1:5]*100, lty=2, lwd=1.5)
legend(x=3.5, y=60, legend=c("Cumulative", "Non-cumulative"),
        lty=c(1,2), pch=c(21,17), bty="n", cex=0.85, pt.bg = "white")
dev.off()

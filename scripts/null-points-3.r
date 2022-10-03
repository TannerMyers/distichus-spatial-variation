#####################################################################################
# This script involves generating null points occurrence points that will be used to
# test if pairs of differentiated lineages of Anolis distichus exhibit niches that 
# are more distinct than expected by chance. Pair 3 is A. d. ignigularis and A. d. properus
# Authors: Tanner C. Myers, Pietro L. H. de Mello, Paul M. Hime, and Richard E. Glor
####################################################################################
# (1) Load packages and assign variables
# Set working directory
working_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation"
setwd(working_dir)

library(tidyverse)
library(data.table)
library(raster)
library(ENMTools)
library(fauxcurrence)

## Command-line args
    args <- commandArgs(trailingOnly = TRUE)
    # seed for reproducibility
    seed <- args[1]

# Load occurrence points for this species pair

# DOUBLE-CHECK THIS PATH
out_dir <- "/scratch/tcm0036/distichus-ddRADseq/analyses/spatial-variation/null-points/" # directory to hold outputs

####################################################################################
# (2) Run fauxcurrence for given seed and species pair

# Prepare data for fauxcurrence
colnames(all_pts) <- c("x", "y", "species")
all_pts$species <- as.character(all_pts$species)

# Use fauxcurrence to obtain bias-corrected background points
set.seed(as.numeric(seed)) # use seed argument as the seed for this run

coords <- as.data.frame()
rast <- raster()
faux_points <- fauxcurrence(coords = coords,
    rast = rast, inter.spp = TRUE,
    sep.inter.spp = TRUE,
    allow.ident.conspec = FALSE)
save(faux_points,
    file = paste0(out_dir, "pair-", pair, "/faux-points-", seed, ".RData"))
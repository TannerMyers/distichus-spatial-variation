library(tidyverse)
library(raster)
library(kuenm)

# Set working directory to directory containing occurence records and input environmental variables
working_dir <- getwd()
setwd(working_dir)

# Load total set of thinned occurrences
occ <- read_csv("locs_joint.csv")

# Load environmental variable PC layers as Raster stack
variables <- raster::stack(list.files(path = "M_variables/Set_1/", pattern = ".asc", full.names = TRUE))

# Define sets of variables to be considered for ENM generation
set <- list(Set_1=c(names(variables)), 
            Set_2=c("pc_1", "pc_2", "pc_3", "pc_4", "pc_5"), 
            Set_3=c("pc_1", "pc_2", "pc_3", "pc_4"))

# `prepare_swd` takes occurrences and raster layers and extracts values for both the occurrence coordinates as well
# as background points that are used as input for Maxent in the SWD format
prepare_swd(occ = occ, species = "lineage", longitude="longitude", latitude="latitude", data.split.method = "random", 
            train.proportion = 0.75, raster.layers=variables, sample.size = 30000, var.sets = set, save = T, 
            name.occ="occurrences", back.folder="M_2", set.seed = 1)

# Assign variables needed for generating candidate models
occ_joint <- "occurrences_joint.csv"
occ_tra <- "occurrences_train.csv"
occ_test <- "occurrences_test.csv"
M_var_dir <- "M_variables"
back_dir <- "M_2"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "no.t.h"
args <- NULL
maxent_path <- "/home/tcm0036/distichus-spatial-variation/enm"
wait <- FALSE
run <- TRUE
kept <- FALSE

#kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
#          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
#          maxent.path = maxent_path, wait = wait, run = run)

out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE

#cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test,
#                        batch = batch_cal, out.eval = out_eval, threshold = threshold,
#                        rand.percent = rand_percent, iterations = iterations, kept = kept,
#                        selection = selection, parallel.proc = paral_proc)

# The function `kuenm_cal_swd` does the entirety of model calibration and evaluation, unlike cases
# in which environmental data is supplied as rasters that executes two functions
kuenm_cal_swd(occ.join = occ_joint, occ.tra = occ_tra, occ.test = occ_test,
             back.dir = back_dir, batch = batch_cal, out.dir.models = out_dir,
             out.dir.eval = out_eval, reg.mult = reg_mult, f.clas = f_clas,
             args = args, maxent.path = maxent_path, selection = selection, kept = kept)
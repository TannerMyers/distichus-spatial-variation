# Set working directory to directory containing occurence records and input environmental variables
working_dir <- "/home/tcm0036/distichus-spatial-variation/enm"
setwd(working_dir)

library(kuenm)

# occ_joint <- "Adist1_all.csv" # All occurrences
# jointPoints <- read.csv(occ_joint)

## load seven environmental variables identified as non-correlated
## CHELSA_bio10_03.asc CHELSA_bio10_04.asc CHELSA_bio10_05.asc CHELSA_bio10_15.asc CHELSA_bio10_16.asc march_EVI_mean.asc may_NDVI_mean.asc
# variables <- raster::stack(list.files(path="M_variables_test/Set_99", pattern=".asc$", full.names=TRUE)) 

# set <- list(Set_1=c(names(variables)), Set_2=c("CHELSA_bio10_03", "CHELSA_bio10_04", "CHELSA_bio10_15", "march_EVI_mean", "may_NDVI_mean"), Set_3=c("CHELSA_bio10_05", "CHELSA_bio10_16","march_EVI_mean", "may_NDVI_mean"),
# 	   Set_4=c("CHELSA_bio10_04","CHELSA_bio10_16","may_NDVI_mean"))

# prepare_swd(occ=jointPoints, species="species", longitude="longitude", latitude="latitude", 
#             data.split.method = "random",
#             train.proportion = 0.5, raster.layers=variables, sample.size = 30000,
#             var.sets = set, save = T, name.occ="occurrences",
#             back.folder="M_2", set.seed = 1)

# Assign variables needed for generating candidate models
occ_joint <- "occurrences_joint.csv"
occ_tra <- "occurrences_train.csv"
occ_test <- "occurrences_test.csv"
back_dir <- "M_2"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
out_dir_eval <- "Calibration_results"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
maxent_path <- "/home/tcm0036/distichus-spatial-variation/enm"
wait <- FALSE
run <- TRUE

## Candidate models are a large set of candidate models created to respond to 
## the need to test broad suites of parameter combinations, such as, distinct regularization multiplier values, 
## various feature classes, and different sets of environmental variables. 
## The following code calls the help page of the function kuenm_cal.
# kuenm_cal_swd(occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
# 	  back.dir = back_dir, batch = batch_cal,
#           out.dir.models = out_dir, out.dir.eval = out_dir_eval,
#           reg.mult = reg_mult, f.clas = f_clas, args = args, max.memory = 1000,
#           maxent.path = maxent_path, selection = "OR_AICc", kept = FALSE)

# Assign variables needed for generating final models
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- FALSE # not projecting model 
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
wait <- FALSE
run <- TRUE
args <- NULL
out.dir <- "mod_swd"

kuenm_mod_swd(occ_joint, back_dir, out_dir_eval, batch_cal, rep_n,
              rep_type, jackknife,
              max.memory = 1000, out_format,
              project, ext_type,
              write_mess, write_clamp, maxent_path,
              args, out.dir, wait, run)

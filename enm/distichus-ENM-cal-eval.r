# Set working directory to directory containing occurence records and input environmental variables
working_dir <- "/home/tcm0036/distichus-spatial-variation/enm"
setwd(working_dir)

library(kuenm)

occ_joint <- "Adist1_all.csv" # All occurrences
jointPoints <- read.csv(occ_joint)

## load seven environmental variables identified as non-correlated
variables1 <- raster::stack(list.files(path="M_variables_test/Set_99", pattern=".asc$", full.names=TRUE))
## load first 5 principal components output by "Raster-PCA.r" script that account for 96.4% of the variation 
variables2 <- raster::stack(list.files(path="M_variables_test/PC_set", pattern=".asc$", full.names=TRUE))
## merge into single stack
variables <- raster::stack(variables1, variables2) 


set <- list(Set_1=c(names(variables1)), Set_2=c(names(variables2)))

prepare_swd(occ=jointPoints, species="species", longitude="longitude", latitude="latitude", 
            data.split.method = "random",
            train.proportion = 0.5, raster.layers=variables2, sample.size = 30000,
            var.sets = set, save = T, name.occ="occurrences",
            back.folder="M_2", set.seed = 1)

# Assign variables
occ_joint <- "occurrences_joint.csv"
occ_tra <- "occurrences_train.csv"
M_var_dir <- "M_variables_test"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
maxent_path <- "/home/tcm0036/distichus-spatial-variation/enm"
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)

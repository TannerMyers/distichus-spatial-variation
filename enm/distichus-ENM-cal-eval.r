# Set working directory to directory containing occurence records and input environmental variables
working_dir <- "/home/tcm0036/distichus-spatial-variation/enm"
setwd(working_dir)

library(kuenm)

#prepare_swd(occ=jointPoints, species="species", longitude="longitude", latitude="latitude", 
#            data.split.method = "random",
#            train.proportion = 0.5, raster.layers=vars, sample.size = 30000,
#            var.sets = set, save = T, name.occ="occurrences",
#            back.folder="M_2", set.seed = 1)

# Assign variables
occ_joint <- "Adist1_all.csv" # All occurrences
occ_tra <- "Adist1_train.csv" # training occurrences to be used for calibration
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

# Set working directory to directory containing occurence records and input environmental variables
working_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/"
setwd(working_dir)

library(kuenm)

# Assign variables
occ_joint <- "Model_calibration/All_Records_30km_thin/Adist1_all.csv" # All occurrences
occ_tra <- "Model_calibration/All_Records_30km_thin/Adist1_train.csv" # training occurrences to be used for calibration
M_var_dir <- "M_variables/"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models/A_distichus_all"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
maxent_path <- working_dir
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
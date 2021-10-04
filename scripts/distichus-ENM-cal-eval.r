
library(kuenm)

# Assign variables
occ_joint <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/Model_calibration/All_Records_30km_thin/Adist1_all.csv" # All occurrences
occ_tra <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/Model_calibration/All_Records_30km_thin/Adist1_train.csv" # training occurrences to be used for calibration
M_var_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/M_variables/"
batch_cal <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/Candidate_models"
out_dir <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/Candidate_Models/A_distichus_all"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
maxent_path <- "/mmfs1/home/tcm0036/distichus-spatial-variation/enm/"
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)
#!/usr/bin/env bash

# Loop over (1) seeds and (2) species pairs to submit jobs to generate null points with "fauxcurrence"
for seed in {1..100}; do
    for pair in 1 2 3 4 5; do
    # pair 1 is A. d. dominicensis 2 and A. d. ignigularis
    # pair 2 is A. d. dominicensis 2 and South paleo-island subspecies
    # pair 3 is A. d. ignigularis and A. d. properus
    # pair 4 is A. d. ignigularis and A. d. ravitergum
    # pair 5 is A. d. ravitergum and South paleo-island subspecies 
            sbatch --job-name==$pair-$seed --partition=bac0071_amd --mem=50G --time=24:00:00 --wrap="

                # Load conda environment with R
                source /mmfs1/home/tcm0036/mambaforge/etc/profile.d/conda.sh 
                conda activate r_env
                
                Rscript null-points-${pair}.r $seed
            "
    done
done
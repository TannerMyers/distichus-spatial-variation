#!/usr/bin/env bash

TAXA=`ls kuenm-analysis`


for species in $TAXA;
	do
		sbatch --job-name=$species-kuenm --partition=bac0071_amd --mail-user=tcm0036@auburn.edu --mail-type=END,FAIL --mem=70G --time=72:00:00 --wrap="
			
			source /mmfs1/home/tcm0036/mambaforge/etc/profile.d/conda.sh 
			conda activate r_env
			
			cd /home/tcm0036/distichus-spatial-variation/enm/kuenm-analysis/$species
					
			/mmfs1/home/tcm0036/mambaforge/envs/r_env/bin/Rscript kuenm-cal-eval.r
		
		"
	done

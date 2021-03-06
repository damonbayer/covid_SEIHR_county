#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 02:00:00   ## 15 min run time limit
#SBATCH --mem=11G 
#SBATCH -o update_projections_4-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //dfs6/pub/bayerd/covid_SEIHR_county

R CMD BATCH --no-save scripts/tidy_pp_gq.R

sbatch --dependency=afterany:$SLURM_JOB_ID slurm_submissions/update_projections_5.sh
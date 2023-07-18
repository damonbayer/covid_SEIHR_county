#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 00:15:00   ## 15 min run time limit
#SBATCH --mem=5G
#SBATCH -o 1_pull_data-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=abakis@uci.edu

module purge
module load R
cd /pub/abakis/GitHub/covid_SEIHR_county

Rscript scripts/pull_case_hospitalizations_data.R

sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/2_fit_model.sh

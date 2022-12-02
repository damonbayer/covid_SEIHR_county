#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 20          ## request 20 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=60G
#SBATCH -o 4_tidy_pp_gq-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //dfs6/pub/bayerd/covid_SEIHR_county

Rscript scripts/tidy_posterior_predictive_and_generated_quantities.R

sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/5a_calcat.sh
sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/5b_figures.sh

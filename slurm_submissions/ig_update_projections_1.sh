#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 00:15:00   ## 15 min run time limit
#SBATCH --mem=5G 
#SBATCH -o update_projections_1-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu

module purge
module load R
cd //pub/igoldst1/covid_SEIHR_county

Rscript scripts/pull_case_hospitalizations_data.R

sbatch --dependency=afterany:$SLURM_JOB_ID ig_update_projections_2.sh

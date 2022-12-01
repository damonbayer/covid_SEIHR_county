#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=3G
#SBATCH -o update_projections_4-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //dfs6/pub/bayerd/covid_SEIHR_county

Rscript scripts/tidy_posterior_predictive_and_generated_quantities.R

# if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
# sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/update_projections_3.sh
# fi

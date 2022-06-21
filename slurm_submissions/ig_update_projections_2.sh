#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=5G 
#SBATCH -o update_projections_2-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=0-56

module purge
module load julia-lts
cd //dfs6/pub/igoldst1/covid_SEIHR_county

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/ig_update_projections_3.sh
fi

julia --project --threads 4 scripts/fit_model_waning.jl $SLURM_ARRAY_TASK_ID
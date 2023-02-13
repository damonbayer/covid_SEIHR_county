#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=5G
#SBATCH -o 2_fit_model-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1-64
#SBATCH --exclude=hpc3-15-29,hpc3-21-30

module purge
module load julia/1.8.5
cd //dfs6/pub/bayerd/covid_SEIHR_county

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/3_gen_pp_gq.sh
fi

julia --project --threads 4 scripts/fit_model.jl $SLURM_ARRAY_TASK_ID
#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 1:00:00   ## 1 hr run time limit
#SBATCH --mem=5G 
#SBATCH -o julia-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1,29,30

module purge
module load julia
cd covid_SEIHR_county

julia --project --threads 4 fit_model.jl $SLURM_ARRAY_TASK_ID

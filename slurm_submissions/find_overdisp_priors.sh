#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=4G
#SBATCH -o update_overdisp_priors-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=abakis@uci.edu
#SBATCH --array=1-64
#SBATCH --error=slurm-out/%x.%A.err
#SBATCH --out=slurm-out/%x.%A.out

module purge
module load R
cd /pub/abakis/git/covid_SEIHR_county

Rscript --no-save  --no-save scripts/find_overdisp_priors.R $SLURM_ARRAY_TASK_ID

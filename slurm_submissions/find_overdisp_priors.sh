#!/bin/bash

#SBATCH -p stats.p
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH --mem=4G
#SBATCH --mail-type=end
#SBATCH --mail-user=abakis@uci.edu
#SBATCH --error=log/%x.%A.err
#SBATCH --out=log/%x.%A.out
#SBATCH --array=1-64

#module purge
#module load R/4.2.1
cd /home/abakis/git/covid_SEIHR_county

Rscript --no-save  --no-save scripts/find_overdisp_priors.R $SLURM_ARRAY_TASK_ID

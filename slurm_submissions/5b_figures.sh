#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=3G
#SBATCH -o 5b_figures-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=abakis@uci.edu
#SBATCH --error=slurm-out/%x.%A.err
#SBATCH --out=slurm-out/%x.%A.out

module purge
module load R
cd /pub/abakis/GitHub/covid_SEIHR_county

Rscript scripts/make_figures.R

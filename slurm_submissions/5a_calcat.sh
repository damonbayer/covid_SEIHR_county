#!/bin/bash

#SBATCH -p stats.p
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mem=3G
#SBATCH --error=log/%x.%A.err
#SBATCH --out=log/%x.%A.out

#module purge
#module load R/4.2.1
cd /home/abakis/git/covid_SEIHR_county

Rscript scripts/format_results_for_calcat.R

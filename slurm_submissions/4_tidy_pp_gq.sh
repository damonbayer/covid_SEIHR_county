#!/bin/bash

#SBATCH -p stats.p
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 04:00:00
#SBATCH --mem=60G
#SBATCH --error=log/%x.%A.err
#SBATCH --out=log/%x.%A.out

#module purge
#module load R/4.2.1
cd /home/abakis/git/covid_SEIHR_county

Rscript scripts/tidy_posterior_predictive_and_generated_quantities.R

sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/5a_calcat.sh
sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/5b_figures.sh

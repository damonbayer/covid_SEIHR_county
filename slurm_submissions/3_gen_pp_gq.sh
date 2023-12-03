#!/bin/bash

#SBATCH -p stats.p
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mem=3G
#SBATCH --error=log/%x.%A.err
#SBATCH --out=log/%x.%A.out
#SBATCH --array=1-64

#module purge
cd /home/abakis/git/covid_SEIHR_county

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/4_tidy_pp_gq.sh
fi

julia +1.8 --project scripts/generate_posterior_predictive_and_generated_quantities.jl $SLURM_ARRAY_TASK_ID

#!/bin/bash

#SBATCH -p stats.p
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 72:00:00
#SBATCH --mem=5G
#SBATCH --mail-type=end
#SBATCH --mail-user=abakis@uci.edu
#SBATCH --error=log/%x.%A.err
#SBATCH --out=log/%x.%A.out
#SBATCH --array=1-64
#SBATCH --exclude=stats-5

#module purge
cd /home/abakis/git/covid_SEIHR_county

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/3_gen_pp_gq.sh
fi

julia +1.8 --project --threads 4 scripts/fit_model.jl $SLURM_ARRAY_TASK_ID

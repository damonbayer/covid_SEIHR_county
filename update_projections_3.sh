#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 00:20:00   ## 15 min run time limit
#SBATCH --mem=11G 
#SBATCH -o julia-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //data/homezvol2/bayerd/covid_SEIHR_county/

R CMD BATCH --no-save make_plots.R

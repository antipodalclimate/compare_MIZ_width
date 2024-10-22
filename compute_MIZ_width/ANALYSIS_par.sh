#!/bin/bash

#SBATCH -t 4:00:00
#SBATCH -n 8
#SBATCH --mem=0
#SBATCH --account=epscor-condo
#SBATCH --mail-user=bndnchrs@gmail.com 
#SBATCH --mail-type=ALL
#SBATCH -J=bybeam
# #SBATCH --constraint=skylake
# #SBATCH --exclusive

matlab-threaded -r drive_AT_stats_par




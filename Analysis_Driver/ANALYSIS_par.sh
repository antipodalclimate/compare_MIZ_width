#!/bin/bash

#SBATCH -t 8:00:00
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem=0
#SBATCH --account=epscor-condo
#SBATCH --mail-user=bndnchrs@gmail.com 
#SBATCH --mail-type=ALL
#SBATCH -J=bybeam
# #SBATCH --constraint=skylake
# #SBATCH --exclusive

module load matlab
matlab-threaded -r drive_MIZ_width


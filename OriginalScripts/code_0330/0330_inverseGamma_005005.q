#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=2
#SBATCH --job-name=InvGama005005

R CMD BATCH ./0330_inverseGamma_005005.R
#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=2
#SBATCH --job-name=InvGama11

R CMD BATCH ./0330_inverseGamma_11.R
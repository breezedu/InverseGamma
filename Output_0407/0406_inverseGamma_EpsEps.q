#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=10
#SBATCH --job-name=InvGama11

R CMD BATCH ./0405_inverseGamma_11.R
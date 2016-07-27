#!/bin/sh

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=16
#SBATCH --job-name=Rstan

R CMD BATCH ./inverse_epseps_with_summary.R
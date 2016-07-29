#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -c 12
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=3
#SBATCH --job-name=lmer4_SQ
R CMD BATCH ./0725_lmer_shuaiqi.R
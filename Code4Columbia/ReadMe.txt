#####################################################
## Read Me 
######################################################


Scripts and dataset in this repository:
Part I: RStan
1. inverse_epseps_with_summary.q
	--- This is the sbatch document, to submit to the cluster.
2. inverse_epseps_with_summary.R 
	--- This is the R script to run the simulation.
3. data/exon_level_process_v2.zip
	--- This is the archived dataset. 
	Before running the R script, make sure to unzip the document.
	After unzipping, a txt document exon_level_process_v2.txt is expected.


####################################################
Pre-required packages for RStan:
R version 3.2.4 or above
Rcpp, RcppEigen, BH, stats4, EDISON
Have to install RStan package before submitting the scripts



#####################################################
Procedure for submitting sbatch jobs to cluster
We use sbatch command to call inverse_epseps_with_summary.q
	gd44@dscr-slogin-01  ~ $ sbatch inverse_epseps_with_summary.q

The sbatch script will allocate memory needed and the R script


######################################################
##	sbatch sample:
## 	Here in this sbatch script, we allocated 16 nodes, 
##      For each node, we allocated 12 CPUs.
## 	For each CPU we could allocate 8G memory.

	#!/bin/sh

	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=jeff.du@duke.edu
	#SBATCH -c 12
	#SBATCH --mem-per-cpu=8G
	#SBATCH --nodes=16
	#SBATCH --job-name=Rstan

	R CMD BATCH ./inverse_epseps_with_summary.R
#####################################################

#####################################################


Part II: LME4

Scripts and dataset in this repository:
4. 0725_lmer_shuaiqi.q
	--- This is the sbatch document, to submit slurm jobs to the cluster.
5. 0725_lmer_shuaiqi.R
	--- This is the R script to run the lmer-simulation.
6. data/exon_level_process_v2.zip
	--- Part II shares the same dataset with Part I, the archieved dataset.
	--- Unzip the file, get exon_level_process_v2.txt document.
	--- 0725_lmer_shuaiqi.R will read data from the txt document.


#####################################################
Pre-required packages for lme4
	--- Linear Mixed-Effects Models using 'Eigen' and S4

R version 3.3 or above (we have tried, it did not work on R 3.2)
Rcpp, knitr, ggplot2, and gamm4. 

#####################################################

#####################################################
Procedure for submitting sbatch jobs to cluster
We use sbatch command to call 0725_lmer_shuaiqi.q
	gd44@dscr-slogin-01  ~ $ sbatch 0725_lmer_shuaiqi.q

The sbatch script will allocate memory needed and the R script


######################################################
##	sbatch sample:
## 	Here in this sbatch script, we allocated 16 nodes, 
##      For each node, we allocated 12 CPUs.
## 	For each CPU we could allocate 8G memory.

	#!/bin/sh

	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=jeff.du@duke.edu
	#SBATCH -c 12
	#SBATCH --mem-per-cpu=8G
	#SBATCH --nodes=5
	#SBATCH --job-name=lme4_SQ

	R CMD BATCH ./0725_lmer_shuaiqi.R
#####################################################

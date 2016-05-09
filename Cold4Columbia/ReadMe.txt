#####################################################
## Read Me 
######################################################


Scripts and dataset in this repository:
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
Rcpp, RcppEigen, BH, stats4, 



#####################################################
Procedure for submitting sbatch jobs to cluster
We use sbatch command to call inverse_epseps_with_summary.q
	gd44@dscr-slogin-01  ~ $ sbatch inverse_epseps_with_summary.q

The sbatch script will allocate memory needed and the R script


######################################################
##	sbatch sample:
## 	Here in this sbatch script, we allocated 16 nodes, 
##      there are 12 cpus on each node, and for each cpu there are 8G memory.

	#!/bin/sh

	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=jeff.du@duke.edu
	#SBATCH -c 12
	#SBATCH --mem-per-cpu=8G
	#SBATCH --nodes=16
	#SBATCH --job-name=Rstan

	R CMD BATCH ./inverse_epseps_with_summary.R
#####################################################

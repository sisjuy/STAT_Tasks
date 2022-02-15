#!/bin/sh
#PBS -l walltime=2400:00:00
cd $PBS_O_WORKDIR
/NA2/R-3.4.4/bin/R CMD BATCH biased_cluster_combine_thetafunc.R 

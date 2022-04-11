#!/bin/bash
#SBATCH --job-name=predict_hybrid_genotype
#SBATCH --nodes=1
#SBATCH --time=0-5:00:00
#SBATCH --account=modelscape
#SBATCH --mem-per-cpu=500G

module load swset/2018.05 gcc/7.3.0 r/4.0.5-py27

R CMD BATCH /gscratch/emcfarl2/predicting_hybrids/likelihood_of_genotypes_WIP_310322.R

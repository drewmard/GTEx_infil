#!/bin/bash -l
#SBATCH -J EBM
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --array=1-189:1

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# Wrapper script for combining related p-values using empirical brown's method

spack load -r r@3.5.0
i=$SLURM_ARRAY_TASK_ID
Rscript /athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/scripts/GeneticAnalysis/shuffled_data/Empirical_Brown_pval_shuffled.R $i

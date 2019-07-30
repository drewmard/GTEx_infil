#!/bin/bash -l
#SBATCH -J EBM
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --array=70-72:1

spack load -r r@3.5.0
i=$SLURM_ARRAY_TASK_ID
Rscript /athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/scripts/GeneticAnalysis/Empirical_Brown_pval.R $i

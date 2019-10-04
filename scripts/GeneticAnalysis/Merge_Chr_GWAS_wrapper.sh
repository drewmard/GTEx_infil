#!/bin/bash -l
#SBATCH -J ChrMg
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --array=1-189:1

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# Script for calling R script that combines GWAS data across chromosomes
# (for cluster computing)

spack load -r r@3.5.0
i=$SLURM_ARRAY_TASK_ID
dir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS
f1=$dir/GTEx.pheno$i.1.qassoc
f2=$dir/GTEx.pheno$i.2.qassoc
f3=$dir/GTEx.pheno$i.3.qassoc

Rscript /athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/scripts/GeneticAnalysis/Merge_Chr_GWAS.R $i

# done
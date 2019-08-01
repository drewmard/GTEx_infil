#!/bin/bash -l
#SBATCH -J GTEx
#SBATCH -n 1
#SBATCH --tasks-per-node=1
#SBATCH --mem=32G
#SBATCH --array=14-17:1

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# Script to run GWAS for identifying iQTLs

source activate HLMM

i=$SLURM_ARRAY_TASK_ID
geno=/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name
pheno=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt
# pheno=/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name.tmp.txt
outdir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS

mkdir -p $outdir

plink --bfile $geno --pheno $pheno --all-pheno --assoc --out $outdir/GTEx.chr$i --allow-no-sex --maf 0.05 --chr $i

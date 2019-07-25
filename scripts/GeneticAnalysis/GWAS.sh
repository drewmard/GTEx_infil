#!/bin/bash -l
#SBATCH -J GTEx
#SBATCH -n 1
#SBATCH --tasks-per-node=4
#SBATCH --mem=4G
#SBATCH --array=1-14:1

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# Script to run GWAS for identifying iQTLs

source activate HLMM

i=$SLURM_ARRAY_TASK_ID
geno=/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name
pheno=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt
outdir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS

mkdir -p $outdir

plink --bfile $geno --pheno $pheno --all-pheno --assoc --out $outdir/GTEx.pheno --allow-no-sex --maf 0.05 --chr $i

# if [ ! -f $outdir/GTEx.pheno$i.1.qassoc ]; then
# plink --bfile $geno --pheno $pheno --pheno-name pheno$i.1 --assoc --out $outdir/GTEx.pheno$i.1 --allow-no-sex --maf 0.05 --chr $i
# fi
# if [ ! -f $outdir/GTEx.pheno$i.2.qassoc ]; then
# plink --bfile $geno --pheno $pheno --pheno-name pheno$i.2 --assoc --out $outdir/GTEx.pheno$i.2 --allow-no-sex --maf 0.05
# fi
# if [ ! -f $outdir/GTEx.pheno$i.3.qassoc ]; then
# plink --bfile $geno --pheno $pheno --pheno-name pheno$i.3 --assoc --out $outdir/GTEx.pheno$i.3 --allow-no-sex --maf 0.05 --chr $i
# fi

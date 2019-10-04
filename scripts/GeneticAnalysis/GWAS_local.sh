#!/bin/bash -l

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# Script to run GWAS for identifying iQTLs

source activate HLMM

# i=$1

for i in {1..22}
do

geno=/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name
pheno=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt
outdir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS

mkdir -p $outdir

plink --bfile $geno --pheno $pheno --all-pheno --assoc --out $outdir/GTEx.chr$i --allow-no-sex --maf 0.05 --chr $i

done

# done
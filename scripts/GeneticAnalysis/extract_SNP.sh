#!/bin/bash -l
#SBATCH -J snp

#source activate GTEx
indir=/athena/elementolab/scratch/anm2868/GTEx
outdir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/Genotype_Subset
prefix=gtex_all.filter.name

mkdir -p $outdir
snp=$1 # 1_751756_T_C_b37
plink=~/miniconda3/envs/GTEx/bin/plink
$plink --bfile $indir/$prefix --snp $snp --recodeA --out $outdir/${prefix}.$snp


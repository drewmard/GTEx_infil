#!/bin/bash -l
#$ -N pca
#$ -j y
#$ -m a
#$ -l h_rt=24:00:00
#$ -pe smp 2
#$ -l h_vmem=4g
#$ -l athena=true

source activate GTEx

TIS=$1

DIR=/athena/elementolab/scratch/anm2868/GTEx
FILE=gtex_all.filter.name.$TIS

# Run PCA
plink --bfile $DIR/$FILE --pca 10 --out $DIR/GENO_PCA/$FILE

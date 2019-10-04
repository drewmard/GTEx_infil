# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu
# script for calculating PCs based on genetic data

source activate GTEx

DIR=/athena/elementolab/scratch/anm2868/GTEx
FILE=gtex_all.filter.name

# Run PCA
plink --bfile $DIR/$FILE --pca 10 --out $DIR/GENO_PCA/$FILE

# done

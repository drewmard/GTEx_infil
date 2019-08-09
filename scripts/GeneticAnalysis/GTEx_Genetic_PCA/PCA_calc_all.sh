source activate GTEx

DIR=/athena/elementolab/scratch/anm2868/GTEx
FILE=gtex_all.filter.name

# Run PCA
plink --bfile $DIR/$FILE --pca 10 --out $DIR/GENO_PCA/$FILE

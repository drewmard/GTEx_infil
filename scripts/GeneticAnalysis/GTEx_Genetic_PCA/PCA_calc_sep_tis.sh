source activate GTEx

DIR=/athena/elementolab/scratch/anm2868/GTEx
FILE=gtex_all.filter.name

# Run PCA

while read -r tis;
do

plink --bfile $DIR/$FILE --pca 10 --keep /athena/elementolab/scratch/anm2868/GTEx/${tis}_ID.txt --out $DIR/GENO_PCA/$FILE.$tis

done < /athena/elementolab/scratch/anm2868/GTEx/tis.uniq.txt


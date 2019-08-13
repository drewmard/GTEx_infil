#!/bin/bash -l

geno=/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name
outdir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS
COLUMN=Pval_Brown

source activate GTEx

for i in {1..189};
do
results=$outdir/GTEx.pheno$i.ALL_EBM.sig.txt
outFile=$outdir/GTEx.pheno$i.ALL_EBM.CLUMPED_1e-5.txt
plink --bfile $geno --clump $results --clump-p1 1e-5 --clump-p2 1e-5 --clump-r2 0.01 --clump-kb 5000 --clump-field $COLUMN --out $outFile
outFile=$outdir/GTEx.pheno$i.ALL_EBM.CLUMPED_5e-8.txt
plink --bfile $geno --clump $results --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 5000 --clump-field $COLUMN --out $outFile
done

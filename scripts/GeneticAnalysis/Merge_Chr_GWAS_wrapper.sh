#!/bin/bash -l
#SBATCH -J ChrMg
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --array=1-221:1

spack load -r r@3.5.0
i=$SLURM_ARRAY_TASK_ID
dir=/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS
f1=$dir/GTEx.pheno$i.1.qassoc
f2=$dir/GTEx.pheno$i.2.qassoc
f3=$dir/GTEx.pheno$i.3.qassoc
if [ ! -f $f1 ] || [ ! -f $f2 ] || [ ! -f $f3 ]; then
Rscript /athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/scripts/GeneticAnalysis/Merge_Chr_GWAS.R $i
else
rm $dir/GTEx.chr*.pheno$i.*.qassoc
fi


# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for combining GWAS data across chromosomes

# Packages
library(data.table)

# Arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) # what pheno to look at?

for (j in 1:3) {
  print(j)
  for (CHR in 1:22) {
    print(CHR)
    df.tmp <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.chr',CHR,'.pheno',i,'.',j,'.qassoc'),data.table=F,stringsAsFactors = F)
    
    if (CHR==1) {
      df.save <- df.tmp
    } else {
      df.save <- rbind(df.save,df.tmp)
    }
  }
  print('Saving...')
  f.save <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.',j,'.qassoc')
  fwrite(df.save,f.save,col.names = T,row.names = F,na='NA',quote=F,sep='\t')
}

# done




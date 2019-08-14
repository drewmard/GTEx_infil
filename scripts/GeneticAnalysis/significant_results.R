# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to extract significant results

library(data.table)
infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
save_every_iteration <- FALSE # TRUE

for (i in 1:nrow(infil_pheno)) {
# for (i in 1:189) {

  print(i)
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.sig.txt')
  if (file.exists(f)) {
    df <- fread(f,data.table = F,stringsAsFactors = F)
    df$tissue <- infil_pheno$tissue[i]
    df$cell <- infil_pheno$cell[i]
    df.sub <- subset(df,Pval_Brown < 5e-8)
    df.sub2 <- subset(df,Pval_Brown < 5e-8 | Pval_Brown_Abs < 5e-8 | p_lrt_cibersort_rel < 5e-8 | p_lrt_cibersort_abs < 5e-8 | p_lrt_xCell < 5e-8 )
    
    if (i == 1) {
      df.sub.save <- df.sub
      df.sub2.save <- df.sub2
    } else {
      df.sub.save <- rbind(df.sub.save,df.sub)
      df.sub2.save <- rbind(df.sub2.save,df.sub2)
    }
    
    if (save_every_iteration) {
      f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig.txt')
      fwrite(df.sub.save,f,col.names = T,row.names = F,sep='\t',quote = F)
      f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig2.txt')
      fwrite(df.sub2.save,f,col.names = T,row.names = F,sep='\t',quote = F)
    }
  } else {
    print(paste0('missing: ',i))
  }
  
}

f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig.txt')
fwrite(df.sub.save,f,col.names = T,row.names = F,sep='\t',quote = F)
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig2.txt')
fwrite(df.sub2.save,f,col.names = T,row.names = F,sep='\t',quote = F)

#########
# analyze in R
library(data.table)
f <- '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig2.txt'
df.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
df.sub.save$phenotype <- paste(df.sub.save$tissue,df.sub.save$cell,sep = '-')
df.sub.save[order(df.sub.save$Pval_Brown,decreasing = F),][1:20,]

x <- df.sub.save
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')
x <- subset(df.sub.save,Pval_Brown < 5e-8)
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')
x <- subset(df.sub.save,Pval_Brown_Abs < 5e-8)
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')
x <- subset(df.sub.save,p_lrt_cibersort_rel < 5e-8)
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')
x <- subset(df.sub.save,p_lrt_cibersort_abs < 5e-8)
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')
x <- subset(df.sub.save,p_lrt_xCell < 5e-8)
paste0(length(unique(paste(x$tissue,x$cell,sep = '-'))),' phenotypes in ',length(unique(x$tissue)), ' tissues.')

# note: output from eQTL_network_gen.R
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.txt')
df.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
df.sub.save$phenotype <- paste(df.sub.save$tissue,df.sub.save$cell,sep = '-')
df.sub.save[order(df.sub.save$Pval_Brown,decreasing = F),][1:20,]


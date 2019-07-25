# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

library(data.table)
infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)

for (i in 1:72) {
  print(i)
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.sig.txt')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  df$tissue <- infil_pheno$tissue[i]
  df$cell <- infil_pheno$cell[i]
  df.sub <- subset(df,Pval_Brown < 5e-8)
  df.sub2 <- subset(df.sub,Pval_Brown < 5e-8 | Pval_Brown_Abs < 5e-8 | p_lrt_cibersort_rel < 5e-8 | p_lrt_cibersort_abs < 5e-8 | p_lrt_xCell < 5e-8 )

  if (i == 1) {
    df.sub.save <- df.sub
    df.sub2.save <- df.sub2
  } else {
    df.sub.save <- rbind(df.sub.save,df.sub)
    df.sub2.save <- rbind(df.sub2.save,df.sub2)
  }
}

f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig.txt')
fwrite(df.sub.save,f,col.names = T,row.names = F,sep='\t',quote = F)
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig2.txt')
fwrite(df.sub2.save,f,col.names = T,row.names = F,sep='\t',quote = F)

f <- '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig.txt'
df.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
df.sub.save[order(df.sub.save$Pval_Brown,decreasing = F),][1:5,]
unique(paste(df.sub.save$tissue,df.sub.save$cell,sep = '-'))
df.min <- aggregate(df.sub.save$Pval_Brown,by=list(tissue=df.sub.save$tissue,cell=df.sub.save$cell),min)
merge(df.min,df.sub.save[,c('tissue','cell','SNP','Pval_Brown')],by.x=c('tissue','cell','x'),by.y=c('tissue','cell','Pval_Brown'))

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to merge clumped results

library(data.table)
infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
save_every_iteration <- FALSE # TRUE
thres <- '5e-8'
for (thres in c('1e-5','5e-8')) {
  FIRST <- TRUE
  for (i in 1:223) {
    print(i)
    f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.CLUMPED_',thres,'.txt.clumped')
    if (file.exists(f)) {
      
      df <- fread(f,data.table = F,stringsAsFactors = F)
      df$tissue <- infil_pheno$tissue[i]
      df$cell <- infil_pheno$cell[i]
      
      if (FIRST) {
        df.save <- df
        FIRST <- FALSE
      } else {
        df.save <- rbind(df.save,df)
      }
      
      # if (save_every_iteration) {
      #   f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.CLUMPED_5e-8.txt.clumped')
      #   fwrite(df.save,f,col.names = T,row.names = F,sep='\t',quote = F)
      # }
    } else {
      print(paste0('missing: ',i))
    }
  }
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.ALL_EBM.CLUMPED_',thres,'.txt.clumped')
  fwrite(df.save,f,col.names = T,row.names = F,sep='\t',quote = F)
}

# phenotypes w/ multiple SNP effects
# df <- fread(f,data.table = F,stringsAsFactors = F)
# sort(table(paste(df$tissue,df$cell)))
# subset(df,tissue=='Lung' & cell=='CD4_Tcells')



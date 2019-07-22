# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for pre-processing of cibersort output

library(data.table)

for (absolute in c('F','T')) {

  # covariate data:
  workdir <- '/athena/elementolab/scratch/anm2868/GTEx/COVAR/'
  df.cov <- fread(paste0(workdir,'GTEx_v7_Annotations_SubjectPhenotypesDS.txt'),data.table = F,stringsAsFactors = F)
  df.tis <- fread(paste0(workdir,'GTEx_v7_Annotations_SampleAttributesDS.txt'),data.table = F,stringsAsFactors = F)
  
  # deconvolution data:
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-',absolute,'.QN-F.perm-1000.DECONV_ONLY.txt')
  df.ciber <- fread(f,data.table = F,stringsAsFactors = F)
  
  # merge together deconvolution data w/ covariate data
  df.ciber.2 <- merge(df.ciber,df.tis[,c('SAMPID','SMTSD')],by.x='Input Sample',by.y = 'SAMPID')
  paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
  df.ciber.2$ID <- sapply(strsplit(df.ciber.2[,1],"-"),paste.s)
  df.ciber.2 <- merge(df.ciber.2, df.cov,by.x = "ID", by.y = 'SUBJID')
  
  # merge immune cell types into broader lineages
  s <- 'CD4_Tcells'
  df.ciber.2[,s] <- apply(df.ciber.2[,grepl('CD4',colnames(df.ciber.2))],1,sum)
  s <- 'MacrophageSum'
  df.ciber.2[,s] <- apply(df.ciber.2[,grepl('Macrophage',colnames(df.ciber.2))],1,sum)
  
  # save
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-',absolute,'.QN-F.perm-1000.txt')
  fwrite(df.ciber.2,f,col.names = T,row.names = F,sep='\t',quote=F,na='NA')
  
}


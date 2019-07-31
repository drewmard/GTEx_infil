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
  s <- 'Bcellsum'
  df.ciber.2[,s] <- apply(df.ciber.2[,c(grep('B cells',colnames(df.ciber.2)),grep('Plasma cells',colnames(df.ciber.2)))],1,sum)
  s <- 'NK_Sum'
  df.ciber.2[,s] <- apply(df.ciber.2[,grepl('NK cells',colnames(df.ciber.2))],1,sum)
  s <- 'DendriticSum'
  df.ciber.2[,s] <- apply(df.ciber.2[,grepl('Dendritic',colnames(df.ciber.2))],1,sum)
  s <- 'MastSum'
  df.ciber.2[,s] <- apply(df.ciber.2[,grepl('Mast cells',colnames(df.ciber.2))],1,sum)
  # s <- 'TcellSum'
  # df.ciber.2[,s] <- apply(df.ciber.2[,grepl('T cells',colnames(df.ciber.2))],1,sum)
  s <- 'Myeloid_Sum'
  df.ciber.2[,s] <- apply(df.ciber.2[,c('MacrophageSum','Neutrophils','DendriticSum','MastSum','Monocytes','Eosinophils')],1,sum)
  #
  s <- 'Lymph_Sum'
  TcellSum <- colnames(df.ciber.2)[grepl('T cells',colnames(df.ciber.2))]
  df.ciber.2[,s] <- apply(df.ciber.2[,c('NK_Sum',TcellSum,'Bcellsum')],1,sum)
  #
  s <- 'CD4.CD8'
  df.ciber.2[,'CD4.CD8'] <- (df.ciber.2[,'CD4_Tcells']+1e-10)/(df.ciber.2[,'T cells CD8']+1e-10)
  df.ciber.2$CD4.CD8[which(df.ciber.2[,'CD4_Tcells']==0 | df.ciber.2[,'T cells CD8']==0)] <- NA
  source('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/rntransform.R')
  for (tis in unique(df.ciber.2$SMTSD)) {
    i <- which(df.ciber.2$SMTSD==tis)
    df.ciber.2$CD4.CD8[i] <- rntransform(df.ciber.2$CD4.CD8[i]) + 10 # +10 such that it passes filter thresholds
  }
  #
  s <- 'Myeloid.Lymph'
  df.ciber.2[,'Myeloid.Lymph'] <- (df.ciber.2[,'Myeloid_Sum']+1e-10)/(df.ciber.2[,'Lymph_Sum']+1e-10)
  df.ciber.2$Myeloid.Lymph[which(df.ciber.2[,'Myeloid_Sum']==0 | df.ciber.2[,'Lymph_Sum']==0)] <- NA
  source('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/rntransform.R')
  for (tis in unique(df.ciber.2$SMTSD)) {
    i <- which(df.ciber.2$SMTSD==tis)
    df.ciber.2$Myeloid.Lymph[i] <- rntransform(df.ciber.2$Myeloid.Lymph[i]) + 10 # +10 such that it passes filter thresholds
  }
  
  # save
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-',absolute,'.QN-F.perm-1000.txt')
  fwrite(df.ciber.2,f,col.names = T,row.names = F,sep='\t',quote=F,na='NA')
  
}


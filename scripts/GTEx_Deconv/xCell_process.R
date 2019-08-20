# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for processing xcell output for downstream analysis

# load
library(data.table)

# covariate data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SubjectPhenotypesDS.txt',data.table = F,stringsAsFactors = F)
df.attr2 <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)
# unique GTEx tissues
tis.uniq = unique(df.attr2$SMTSD)

# Combine different xCell tissue outputs
first <- TRUE
for (tis in tis.uniq) {
  if (first) {
    df.xcell <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out/XCell_',tis,'.txt'),data.table=F)
    df.xcell$SMTSD <- tis
    first <- FALSE
  } else {
    df.tmp <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out/XCell_',tis,'.txt'),data.table=F)
    df.tmp$SMTSD <- tis
    df.xcell <- rbind(df.xcell,df.tmp)
  }
}

# Merge w/ some covariate data
paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
IID <- sapply(strsplit(df.xcell$SAMP,"-"),paste.s)
df.xcell$IID <- IID
df.xcell <- merge(df.xcell,df.attr,by.x='IID',by.y='SUBJID')

# merge related cell types into broader lineages
s <- 'CD4_memory'
df.xcell[,s] <- apply(df.xcell[,c('CD4+ Tcm',
                               'CD4+ Tem',
                               'CD4+ memory T-cells')],1,sum)
s <- 'MacrophageSum'
df.xcell[,s] <- apply(df.xcell[,grepl('Macrophage',colnames(df.xcell))],1,sum)
s <- 'CD8Sum'
df.xcell[,s] <- apply(df.xcell[,grepl('CD8',colnames(df.xcell))],1,sum)
s <- 'Bcellsum'
df.xcell[,s] <- apply(df.xcell[,c(grep('B-cell',colnames(df.xcell)),grep('Plasma cells',colnames(df.xcell)))],1,sum)
s <- 'DendriticSum'
df.xcell[,s] <- apply(df.xcell[,grepl('DC',colnames(df.xcell))],1,sum)
s <- 'Th_Sum'
df.xcell[,s] <- apply(df.xcell[,grepl('Th',colnames(df.xcell))],1,sum)
#
s <- 'Myeloid_Sum'
df.xcell[,s] <- apply(df.xcell[,c('MacrophageSum','Neutrophils','DendriticSum','Mast cells','Monocytes','Eosinophils','Basophils')],1,sum)
#
s <- 'Lymph_Sum'
TcellSum <- colnames(df.xcell)[c(7:15,51,62:65)]
df.xcell[,s] <- apply(df.xcell[,c('NK cells',TcellSum,'Bcellsum')],1,sum)
#
s <- 'CD4Sum'
df.xcell[,s] <- apply(df.xcell[,c('CD4+ Tcm',
                                  'CD4+ Tem',
                                  'CD4+ memory T-cells',
                                  'CD4+ naive T-cells',
                                  'CD4+ T-cells',
                                  'Tregs',
                                  'Th1 cells',
                                  'Th2 cells',
                                  'Tgd cells')],1,sum)
#
s <- 'CD4.CD8'
df.xcell[,'CD4.CD8'] <- (df.xcell[,'CD4Sum']+1e-10)/(df.xcell[,'CD8Sum']+1e-10)
df.xcell$CD4.CD8[which(df.xcell[,'CD4Sum']==0 | df.xcell[,'CD8Sum']==0)] <- NA
source('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/rntransform.R')
for (tis in unique(df.xcell$SMTSD)) {
  i <- which(df.xcell$SMTSD==tis)
  df.xcell$CD4.CD8[i] <- rntransform(df.xcell$CD4.CD8[i]) + 10 # +10 such that it passes filter thresholds
}
#
s <- 'Myeloid.Lymph'
df.xcell[,'Myeloid.Lymph'] <- (df.xcell[,'Myeloid_Sum']+1e-10)/(df.xcell[,'Lymph_Sum']+1e-10)
df.xcell$Myeloid.Lymph[which(df.xcell[,'Myeloid_Sum']==0 | df.xcell[,'Lymph_Sum']==0)] <- NA
source('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/rntransform.R')
for (tis in unique(df.xcell$SMTSD)) {
  i <- which(df.xcell$SMTSD==tis)
  df.xcell$Myeloid.Lymph[i] <- rntransform(df.xcell$Myeloid.Lymph[i]) + 10 # +10 such that it passes filter thresholds
}



# cellTypes <- c('CD8Sum','CD4Sum','Neutrophils','MacrophageSum',
#                'Bcellsum','NK cells','DendriticSum','Mast cells','TcellSum',
#                'Plasma cells','TFH_Sum','Tregs','Tgd cells',
#                'Monocytes','Eosinophils')

# save
fwrite(df.xcell,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',row.names = F,col.names = T,sep = '\t',na='NA',quote=F)



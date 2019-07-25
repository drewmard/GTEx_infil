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
s <- 'CD4Sum'
df.xcell[,s] <- apply(df.xcell[,grepl('CD4',colnames(df.xcell))],1,sum)
s <- 'MacrophageSum'
df.xcell[,s] <- apply(df.xcell[,grepl('Macrophage',colnames(df.xcell))],1,sum)
s <- 'CD8Sum'
df.xcell[,s] <- apply(df.xcell[,grepl('CD8',colnames(df.xcell))],1,sum)
s <- 'Bcellsum'
df.xcell[,s] <- apply(df.xcell[,c(grep('B-cell',colnames(df.xcell)),grep('Plasma cells',colnames(df.xcell)))],1,sum)
s <- 'DendriticSum'
df.xcell[,s] <- apply(df.xcell[,grepl('DC',colnames(df.xcell))],1,sum)
s <- 'TcellSum'
df.xcell[,s] <- apply(df.xcell[,c(7:15,51,62:65)],1,sum)
s <- 'Th_Sum'
df.xcell[,s] <- apply(df.xcell[,grepl('Th',colnames(df.xcell))],1,sum)
s <- 'Lymph_Sum'
df.xcell[,s] <- apply(df.xcell[,c('NK cells','TcellSum','Bcellsum')],1,sum)


# cellTypes <- c('CD8Sum','CD4Sum','Neutrophils','MacrophageSum',
#                'Bcellsum','NK cells','DendriticSum','Mast cells','TcellSum',
#                'Plasma cells','TFH_Sum','Tregs','Tgd cells',
#                'Monocytes','Eosinophils')

# save
fwrite(df.xcell,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',row.names = F,col.names = T,sep = '\t',na='NA',quote=F)



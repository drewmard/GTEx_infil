# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for processing xcell output for downstream analysis

# load
library(data.table)

# covariate data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SubjectPhenotypesDS.txt',data.table = F,stringsAsFactors = F)
# df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)
# unique GTEx tissues
tis.uniq = unique(df.attr$SMTSD)

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
df.xcell$CD4Sum <- (apply(df.xcell[,7:11],1,sum))
df.xcell$MacrophageSum <- (apply(df.xcell[,33:35],1,sum))
df.xcell$CD8Sum <- (apply(df.xcell[,12:15],1,sum))

# save
fwrite(df.xcell,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',row.names = F,col.names = T,sep = '\t',na='NA',quote=F)



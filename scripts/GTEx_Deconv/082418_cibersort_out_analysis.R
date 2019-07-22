library(data.table)
#library(ggplot2)
#library(stringr)
#library(ggpubr)
#library(cowplot)
#library(data.table)

abs <- 'F'
# workdir <- '/Volumes/SeagateBackupPlusDrive/Elemento/'
workdir <- '/athena/elementolab/scratch/anm2868/GTEx/'
df.cov <- fread(paste0(workdir,'GTEx_v7_Annotations_SubjectPhenotypesDS.txt'),data.table = F)
df.tis <- fread(paste0(workdir,'GTEx_v7_Annotations_SampleAttributesDS.txt'),data.table = F)
df.id <- fread(paste0(workdir,'GTEx_v7_genexpr_SPLIT.txt'),data.table = F)

pheno_num <- 1
f <- paste0(workdir,'CIBERSORT_out/','GTEx_v7_genexpr_',pheno_num,'.CIBERSORT.ABS-',abs,'.QN-F.perm-1000.txt')
df.ciber <- fread(f,data.table = F)
#df <- fread('/Volumes/SeagateBackupPlusDrive/Elemento/CIBERSORT.GTEx_v7_genexpr_1.txt',data.table=F)
for (pheno_num in 2:10) {
  f <- paste0(workdir,'CIBERSORT_out/','GTEx_v7_genexpr_',pheno_num,'.CIBERSORT.ABS-',abs,'.QN-F.perm-1000.txt')
  df.ciber.tmp <- fread(f,data.table = F)
  df.ciber <- rbind(df.ciber,df.ciber.tmp)
}

df.ciber.2 <- merge(df.ciber,df.tis[,c('SAMPID','SMTSD')],by.x='Input Sample',by.y = 'SAMPID')
paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
df.ciber.2$ID <- sapply(strsplit(df.ciber.2[,1],"-"),paste.s)
df.ciber.2 <- merge(df.ciber.2, df.cov,by.x = "ID", by.y = 'SUBJID')
df.ciber.2$AGE <- as.factor(df.ciber.2$AGE)

# CD4
s <- 'CD4_Tcells'
df.ciber.2[,s] <- apply(df.ciber.2[,grepl('CD4',colnames(df.ciber.2))],1,sum)
s <- 'MacrophageSum'
df.ciber.2[,s] <- apply(df.ciber.2[,grepl('Macrophage',colnames(df.ciber.2))],1,sum)

f <- paste0(workdir,'CIBERSORT_out/','GTEx_v7_genexpr_','ALL','.CIBERSORT.ABS-',abs,'.QN-F.perm-1000.txt')
fwrite(df.ciber.2,f)




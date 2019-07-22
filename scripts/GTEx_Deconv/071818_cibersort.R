# Install
#install.packages('e1071',lib='/home/anm2868/lib',repos='https://mirrors.sorengard.com/cran/')
#source("http://bioconductor.org/biocLite.R")
#biocLite("preprocessCore",lib='/home/anm2868/lib')
#install.packages('data.table',lib='/home/anm2868/lib',repos='https://mirrors.sorengard.com/cran/')

# load 
library('e1071') #,lib.loc='/home/anm2868/lib')
library('parallel')
library('preprocessCore') #,lib.loc='/home/anm2868/lib')
library('data.table') #,lib.loc='/home/anm2868/lib')

args = commandArgs(trailingOnly=TRUE)
absolute=args[1]
print(absolute)
#split_num=args[2]
#print(split_num)

library('data.table')
source('/home/anm2868/scripts/CIBERSORT.R')
sig.matrix.file<-"/home/anm2868/etc/LM22.txt"
#for (i in split_num:split_num) {
for (i in 5:10) {
#for (i in 1:1) {
  print(paste0('Running ',i,'...'))
  
  mixture.file<-paste0("/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_v7_genexpr_",i,".txt")
  
  if (absolute=='abs') {
    results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=TRUE, abs_method='sig.score')
  } else {
    results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=FALSE)
  }
  results <- as.data.frame(results)
  results$`Input Sample` <- rownames(results)
  results <- results[,c(ncol(results),1:(ncol(results)-1))]
  rownames(results) <- NULL
  
  if (absolute=='abs') {
    fwrite(results,file=paste0("/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/GTEx_v7_genexpr_",i,".CIBERSORT.ABS-T.QN-F.perm-1000.txt"),quote=F,row.names=F,sep='\t')
  } else {
    fwrite(results,file=paste0("/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/GTEx_v7_genexpr_",i,".CIBERSORT.ABS-F.QN-F.perm-1000.txt"),quote=F,row.names=F,sep='\t')
  }
}

# # For RELATIVE p-values
# results <- CIBERSORT(sig.matrix.file,mixture.file, perm=0, QN=FALSE, absolute=FALSE)

# example output from CIBERSORT server
# X <- fread('/Volumes/SeagateBackupPlusDrive/Elemento/CIBERSORT.GTEx_v7_genexpr_1.txt')

# example R output from CIBERSORT server
# X <- fread(paste0("/Volumes/SeagateBackupPlusDrive/Elemento/CIBERSORT_out/GTEx_v7_genexpr_",i,".CIBERSORT.ABS-T.QN-F.perm-0.txt"),sep='\t')

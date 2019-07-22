# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

library('e1071') #,lib.loc='/home/anm2868/lib')
library('parallel')
library('preprocessCore') #,lib.loc='/home/anm2868/lib')
library('data.table') #,lib.loc='/home/anm2868/lib')

# Inputs
args <- c('abs')
# args = commandArgs(trailingOnly=TRUE)
absolute=args[1]
print(absolute)
source('/home/anm2868/scripts/CIBERSORT.R')
sig.matrix.file<-"/home/anm2868/etc/LM22.txt"
mixture.file <- '/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/TPM_GeneExpression_for_Mixtures.txt'  
# run cibersort
if (absolute=='abs') {
    results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=TRUE)
    # results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')
} else {
    results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=FALSE)
}
results <- as.data.frame(results)
results$`Input Sample` <- rownames(results)
results <- results[,c(ncol(results),1:(ncol(results)-1))]
rownames(results) <- NULL
  
if (absolute=='abs') {
#    fwrite(results,file=paste0("/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/SynMix_",i,".CIBERSORT.ABS-T.QN-F.perm-1000.txt"),quote=F,row.names=F,sep='\t')
    fwrite(results,file=paste0("/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/SynMix_",i,".CIBERSORT.ABS-T-nosumto1.QN-F.perm-1000.txt"),quote=F,row.names=F,sep='\t')
  } else {
    fwrite(results,file=paste0("/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/SynMix_",i,".CIBERSORT.ABS-F.QN-F.perm-1000.txt"),quote=F,row.names=F,sep='\t')
  }
}


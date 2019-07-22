# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# load ggplot R libraries
library('e1071')
library('parallel')
library('preprocessCore') 
library('data.table') 

# directories: requires path to CIBERSORT R script. need to request from https://cibersort.stanford.edu
CIBERSORT_script <- '/home/anm2868/scripts/CIBERSORT.R'
sig.matrix.file<-"/home/anm2868/etc/LM22.txt"
mixture.file <- '/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'

# command line argument: cibersort rel or abs?
args = commandArgs(trailingOnly=TRUE)
absolute=args[1]
absolute='abs'
# i=args[2]
# print(paste('Running...',absolute,i))
print(paste('Running...',absolute))
source(CIBERSORT_script)

# Run CIBERSORT...
# mixture.file<-paste0("/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_v7_genexpr_",i,".txt")
if (absolute=='abs') {
  results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=TRUE, abs_method='sig.score')
} else if (absolute=='rel') {
  results <- CIBERSORT(sig.matrix.file,mixture.file, perm=1000, QN=FALSE, absolute=FALSE)
}
# format the output:
results <- as.data.frame(results)
results$`Input Sample` <- rownames(results)
results <- results[,c(ncol(results),1:(ncol(results)-1))]
rownames(results) <- NULL

# Save output
if (absolute=='abs') {
  dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles')
  fwrite(results,file='/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.DECONV_ONLY.txt',quote=F,row.names=F,sep='\t',col.names = T)
  
  # create relative CIBERSORT proportions from absolute scores
  results.rel <- results
  for (i in 1:nrow(results)) {
    results.rel[i,2:23] <- results.rel[i,2:23]/results.rel[i,"Absolute score (sig.score)"]
  }
  fwrite(results.rel,file='/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.DECONV_ONLY.txt',quote=F,row.names=F,sep='\t',col.names = T)
  
} else if (absolute=='rel') {
  fwrite(results,file='/athena/elementolab/scratch/anm2868/GTEx/CIBERSORT_out/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.DECONV_ONLY.txt',quote=F,row.names=F,sep='\t',col.names = T)
}




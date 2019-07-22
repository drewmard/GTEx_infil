library(xCell)
library(data.table)
f <- '/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'
print("loading expression...")
df.exp <- fread(f,data.table = F,stringsAsFactors = F)
print("done loading.")
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)

TPM_df.full <- df.exp[,c(2:ncol(df.exp))]
colnames(TPM_df.full)[1] <- 'Gene'
TPM_df.full <- as.data.table(TPM_df.full)
TPM_df.full <- TPM_df.full[ , lapply(.SD, mean), by = list(TPM_df.full$Gene)]
TPM_df.full <- as.data.frame(TPM_df.full)
rownames(TPM_df.full) <- TPM_df.full[,1]
TPM_df.full <- TPM_df.full[,-1]

tis.uniq = unique(df.attr$SMTSD)

##################################################
##################################################
##################################################
# for doing full thing at once
# Result_CellTypes_xCell <- xCellAnalysis(TPM_df.full,rnaseq = T)
# Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
# Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
# 
# fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/infil_output/XCell_','Marderstein','.txt'),sep='\t',quote=F,row.names=F,col.names=T)
##################################################
##################################################
##################################################

##################################################
##################################################
##################################################
# for subsetting by tissue
for (i in 29:length(tis.uniq)) {
  tis = tis.uniq[i]
  print(tis)
  df.sub <- subset(df.attr,SMTSD==tis)
  ind <- which(colnames(TPM_df.full) %in% df.sub$SAMPID)
  if (length(ind)==0) {print(paste0('None: ',tis)); next}
  TPM_df <- TPM_df.full[,ind]

  # Result_CellTypes_xCell <- xCellAnalysis(TPM_df,rnaseq = T)
  Result_CellTypes_xCell <- xCellAnalysis(TPM_df,rnaseq = T,parallel.sz=16)
  Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
  Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
  
  fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/infil_output/xcell_sub_out/XCell_',tis,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
}
##################################################
##################################################
##################################################



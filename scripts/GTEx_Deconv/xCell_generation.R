library(xCell)
library(data.table)
print("loading expression...")
f <- '/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'
df.exp <- fread(f,data.table = F,stringsAsFactors = F)
print("done loading...")

# covariate data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)

# arrange TPM file into right format (merging duplicates into means):
TPM_df.full <- df.exp[,-1] # remove other gene name column
# merge duplicated genes
colnames(TPM_df.full)[1] <- 'Gene'
TPM_df.full <- as.data.table(TPM_df.full); TPM_df.full <- TPM_df.full[ , lapply(.SD, mean), by = list(TPM_df.full$Gene)]; TPM_df.full <- as.data.frame(TPM_df.full)
# and assign gene names to row names
rownames(TPM_df.full) <- TPM_df.full[,1]
TPM_df.full <- TPM_df.full[,-1]

# unique GTEx tissues
tis.uniq = unique(df.attr$SMTSD)

##################################################
##################################################
##################################################
# for performing full xCell at once: NOT RECOMMENDED. PERFORMS WORSE.
# Result_CellTypes_xCell <- xCellAnalysis(TPM_df.full,rnaseq = T)
# Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
# Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
# dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out/xcell_together')
# fwrite(Result_CellTypes_xCell.t,file = '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out/xcell_together/xcell_together.txt',sep='\t',quote=F,row.names=F,col.names=T)
##################################################
##################################################
##################################################

dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out')
# for subsetting by tissue
for (i in 1:length(tis.uniq)) {
  tis = tis.uniq[i]
  print(tis)
  
  # subset:
  df.sub <- subset(df.attr,SMTSD==tis)
  ind <- which(colnames(TPM_df.full) %in% df.sub$SAMPID)
  if (length(ind)==0) {print(paste0('None: ',tis)); next} # debug tool
  TPM_df <- TPM_df.full[,ind]

  # RUN XCELL:
  # Result_CellTypes_xCell <- xCellAnalysis(TPM_df,rnaseq = T)
  Result_CellTypes_xCell <- xCellAnalysis(TPM_df,rnaseq = T,parallel.sz=8) # parallelize
  Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
  Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
  
  # SAVE:
  fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/xcell_sub_out/XCell_',tis,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
}



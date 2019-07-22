# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

library(xCell)
library(data.table)

# Read in synthetic mix TPMs and arrange in format for xCell
f <- '/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/TPM_GeneExpression_for_Mixtures.txt'
TPM_df.full <- fread(f,data.table = F,stringsAsFactors = F)
colnames(TPM_df.full)[1] <- 'Gene'
TPM_df.full <- as.data.table(TPM_df.full)
TPM_df.full <- TPM_df.full[ , lapply(.SD, mean), by = list(TPM_df.full$Gene)]
TPM_df.full <- as.data.frame(TPM_df.full)
rownames(TPM_df.full) <- TPM_df.full[,1]
TPM_df.full <- TPM_df.full[,-1]

##################################################
##################################################
##################################################
# for doing full thing at once. performs worse.
Result_CellTypes_xCell <- xCellAnalysis(TPM_df.full,rnaseq = T)
Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL

fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','','.txt'),sep='\t',quote=F,row.names=F,col.names=T)
##################################################
##################################################
##################################################
# for doing 1 tissue type at a time:
# Tissue 1/2
ind <- which(colnames(TPM_df.full) %in% paste0("Mixture_",1:40,".tab"))
TPM_df.sub <- TPM_df.full[,ind]
Result_CellTypes_xCell <- xCellAnalysis(TPM_df.sub,rnaseq = T)
Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','.sub1','.txt'),sep='\t',quote=F,row.names=F,col.names=T)

# Tissue 2/2
ind <- which(colnames(TPM_df.full) %in% paste0("Mixture_",41:80,".tab"))
TPM_df.sub <- TPM_df.full[,ind]
Result_CellTypes_xCell <- xCellAnalysis(TPM_df.sub,rnaseq = T)
Result_CellTypes_xCell.t <- as.data.frame(t(Result_CellTypes_xCell))
Result_CellTypes_xCell.t$SAMP <- rownames(Result_CellTypes_xCell.t); rownames(Result_CellTypes_xCell.t) <- NULL
fwrite(Result_CellTypes_xCell.t,file = paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','.sub2','.txt'),sep='\t',quote=F,row.names=F,col.names=T)

##################################################
##################################################
##################################################

library(data.table)
# 1 # xcell performed together
# xcell <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','','.txt'),data.table=F,stringsAsFactors = F)

# 2 # xcell performed on tissue subsets
xcell1 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','.sub1','.txt'),data.table=F,stringsAsFactors = F)
xcell2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput','.sub2','.txt'),data.table=F,stringsAsFactors = F)
xcell = rbind(xcell1,xcell2)

# mappings for sample descriptors
mapping <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/etc/SampleSorted_renaming.txt',data.table = F,stringsAsFactors = F,header = F); colnames(mapping)[2] <- 'Long_Name'
synmix <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/etc/mix_prop.txt',data.table = F,stringsAsFactors = F,header=F)
xcell$`SAMP` <- substring(xcell$`SAMP`,1,nchar(xcell$`SAMP`)-4)
xcell <- merge(xcell,mapping,by.x='SAMP',by.y='V1')

# pull out name
x <- strsplit(substring(xcell$Long_Name,1,nchar(xcell$Long_Name)-11),'pct.')

# script to merge cell types
quantify = function(y,code='CD4 Naive') {
  tot = 0
  y[length(y)] <- substring(y[length(y)],1,nchar(y[length(y)])-3)
  for (i in 1:length(y)) {
    y2 <- strsplit(y[i],'_')[[1]]
    celltype <- synmix[synmix[,1]==y2[1],2]
    if (length(celltype) == 0) {next}
    if (celltype==code) {
      tot <- tot + as.numeric(y2[2])
    }
  }
  return(tot)
}

xcell$CD4_Naive_Actual <- as.numeric(lapply(x,quantify,code='CD4 Naive'))
xcell$CD8_Naive_Actual <- as.numeric(lapply(x,quantify,code='CD8 Naive'))
xcell$Tot_Infil_Actual <- xcell$CD4_Naive_Actual + xcell$CD8_Naive_Actual
xcell$Heterogeneity <- as.logical(lapply(x,function(y) length(y) > 3))
xcell$CD4_Tcells_Sum <- apply(xcell[,7:11],1,sum)
xcell$CD8_Tcells_Sum <- apply(xcell[,12:15],1,sum)
xcell$CD4_Naive_Actual_Pct <- xcell$CD4_Naive_Actual/xcell$Tot_Infil_Actual
xcell$CD8_Naive_Actual_Pct <- xcell$CD8_Naive_Actual/xcell$Tot_Infil_Actual

# Save xcell output
fwrite(xcell,'/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput_merged.txt',quote = F,row.names = F,col.names = T,sep = '\t')

# correlations
cor(xcell$CD4_Naive_Actual,xcell$CD4_Tcells_Sum)
cor(xcell$CD4_Naive_Actual_Pct,xcell$CD4_Tcells_Sum)
cor(xcell$CD8_Naive_Actual,xcell$CD8_Tcells_Sum)
cor(xcell$CD8_Naive_Actual_Pct,xcell$CD8_Tcells_Sum)





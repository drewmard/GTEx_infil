# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# What to analyze: 
# based on whether tissue has noticable infiltration across samples, 
# whether cell type is well-represented in samples, 
# and whether sample size is sufficient.

library(data.table)
# df.rel <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table=F,stringsAsFactors = F)
df.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)




tis.uniq <- data.frame(table(df.rel$SMTSD))
df.results <- data.frame(tissue=c(),cell=c())
for (i in 1:length(tis.uniq)) {
  TISSUE <- tis.uniq[i,1]
  df.sub <- subset(df.rel,SMTSD==TISSUE)
  
  # cond 1: does this tissue have sufficient infiltration?
  ciberSigSig <- subset(df.sub, `P-value` < 0.5)
  infil <- (nrow(ciberSigSig)/nrow(df.sub))
  
  # cond 2: is this cell type represented in the sample?
  cellTypes <- c('T cells CD8','CD4_Tcells','Neutrophils','MacrophageSum','T cells regulatory (Tregs)')
  Mc <- as.matrix(ciberSigSig[, cellTypes])
  rownames(Mc) <- ciberSigSig$ID
  cellTypeFreq <- apply(Mc, 2, mean)
  High.Rel <- cellTypes[as.numeric(which(cellTypeFreq > 0.05))]
  
  # cond 3: are there enough samples total?
  N <- nrow(df.sub) # not used, but what was used for GWAS
  
  if (infil > 0.5 & N >= 70) {
    for (cell in High.Rel) {
      df.results <- rbind(df.results,data.frame(tissue=TISSUE,cell=cell))
    }
  }
  
}
df.results$tissue <- as.character(df.results$tissue)
df.results$cell <- as.character(df.results$cell)
fwrite(df.results,'/Volumes/SeagateBackupPlusDrive/Elemento/tissue_x_celltype.txt',sep='\t',quote=F)

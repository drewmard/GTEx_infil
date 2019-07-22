# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# What to analyze: 
# based on whether tissue has noticable infiltration across samples, 
# whether cell type is well-represented in samples, 
# and whether sample size is sufficient.
# ALSO: xcell and cibersort need to "generally" agree w/ each other (no inverse correlation)

library(data.table)
df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)

genetic_data_fam <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/gtex_all.filter.name.fam',data.table = F,stringsAsFactors = F)

# filter: indiv needs to have genetic data
df.rel <- subset(df.rel,ID %in% genetic_data_fam[,2])

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
  N <- tis.uniq[i,2]
  
  if (infil > 0.5 & N >= 70) {
    for (cell in High.Rel) {
      
      # cond 4: do xcell and cibersort "somewhat" agree?
      df.abs.sub <- subset(df.abs,SMTSD==TISSUE)
      df.xcell.sub <- subset(df.xcell,SMTSD==TISSUE)
      
      df.abs.sub[,cell]
      
      df.results <- rbind(df.results,data.frame(tissue=TISSUE,cell=cell))
    }
  }
  
}
df.results$tissue <- as.character(df.results$tissue)
df.results$cell <- as.character(df.results$cell)
fwrite(df.results,'/Volumes/SeagateBackupPlusDrive/Elemento/infiltration_phenotypes.txt',sep='\t',quote=F)

library('data.table')

dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv',showWarnings = F)
dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/Expr_PCA',showWarnings = F)

print('Uploading expression file...')
f <- '/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct'
df.exp <- fread(f,data.table = F,stringsAsFactors = F)
df.tis <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F)

i=0
for (tis in unique(df.tis$SMTSD)) {

  i=i+1
  print(paste0('Beginning tissue ',i,': ',tis))
  
  tis2 <- gsub(' ','_',tis)

  df.sub <- subset(df.tis,SMTSD==tis)
  ind <- which(colnames(df.exp) %in% df.sub$SAMPID)
  df.exp.sub <- df.exp[,c(2,ind)]
  colnames(df.exp.sub)[1] <- 'Gene'
  df.exp.sub <- as.data.table(df.exp.sub); 
  df.exp.sub <- df.exp.sub[ , lapply(.SD, mean), by = list(df.exp.sub$Gene)]; 
  df.exp.sub <- as.data.frame(df.exp.sub)
  rownames(df.exp.sub) <- df.exp.sub[,1]; df.exp.sub <- df.exp.sub[,-1]

  df.exp.sub <- as.data.frame(t(df.exp.sub))

  pca_results <- prcomp(df.exp.sub)
  
  eigs <- pca_results$sdev^2
  perc.var.exp <- eigs[1:4]/sum(eigs)
  perc.var.exp.df <- data.frame(Tissue=tis,PC=paste0('PC',1:4),Percent_Variance_Explained=perc.var.exp)
  if (i==1) {
    perc.var.exp.df.save <- perc.var.exp.df
  } else {
    perc.var.exp.df.save <- rbind(perc.var.exp.df.save, perc.var.exp.df)
  }

  df.PEER <- as.data.frame(pca_results$x[,1:10])
  df.PEER$SAMP <- rownames(df.PEER)
  rownames(df.PEER) <- NULL
  df.PEER <- df.PEER[,c(ncol(df.PEER),1:10)]
  # head(df.PEER)
  # fwrite(df.PEER,paste0('/athena/elementolab/scratch/anm2868/GTEx/PCA/',tis2,'.PCA.txt'),quote = F,row.names = F,sep='\t',col.names = T)
  fwrite(df.PEER,paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/Expr_PCA/',tis2,'.PCA.txt'),quote = F,row.names = F,sep='\t',col.names = T)
}

fwrite(perc.var.exp.df.save,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/Expr_PCA/perc.var.exp.df.save.txt',quote = F,row.names = F,sep='\t',col.names = T)


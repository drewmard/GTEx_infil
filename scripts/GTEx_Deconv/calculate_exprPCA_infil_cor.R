# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: compare PCs of gene expression with deconvolution estimates in each tissue


library('data.table')

infiltration_phenotypes <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
infiltration_phenotypes$tissue2 <- gsub(' ','_',infiltration_phenotypes$tissue)

# infiltration
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

NEW=T
for (i in 1:nrow(infiltration_phenotypes)) {
  s <- infiltration_phenotypes$tissue[i]; s2 <- infiltration_phenotypes$tissue2[i]
  print(paste0('Beginning tissue ',i,': ',s))
  df.PCA = fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/Expr_PCA/',s2,'.PCA.txt'),stringsAsFactors = F,data.table = F)
  df.infil.sub = subset(df.abs,SMTSD==s)
  df.mg <- merge(df.PCA,df.infil.sub,by.x='SAMP',by.y='Input Sample')
  
  cell <- infiltration_phenotypes$cell[i]
  for (PCnum in 2:11) {
    df.tmp = data.frame(tissue=s,cell=cell,PCnum=PCnum-1,pval=summary(lm(df.mg[,PCnum]~df.mg[,cell]))$coef[2,4])
    if (NEW) {
      df.save = df.tmp
      NEW=F
    } else {df.save <- rbind(df.save,df.tmp)}
  }
}


df.sub = subset(df.save,PCnum%in%1:4)
df.sub$pval.fdr = p.adjust(df.sub$pval,method = 'fdr')

df.aggre = aggregate(df.sub[,c('pval','pval.fdr')],by=list(tissue=df.sub$tissue,cell=df.sub$cell),min)
head(df.aggre[order(df.aggre$pval.fdr),])

sum(df.aggre$pval.fdr < 0.1)
nrow(df.aggre)
sum(df.aggre$pval < 0.05)

fwrite(df.aggre,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/PCA_infil_cor.txt',sep='\t',row.names = F,col.names = T,quote=F)
fwrite(df.sub,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GTEx_Deconv/PCA_infil_cor_full.txt',sep='\t',row.names = F,col.names = T,quote=F)


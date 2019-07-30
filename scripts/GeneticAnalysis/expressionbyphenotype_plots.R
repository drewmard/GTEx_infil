# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to draw genotype by phenotype plots



library(data.table)
library(stringr)
library(ggplot2)

# df.geno <- fread('/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',data.table = F,stringsAsFactors = F)

tis <- 'Whole Blood'

# tis <- 'Colon - Sigmoid'
# GENE <- 'C22orf43'
# CELL <- 'Lymph_Sum'

# tis <- 'Heart - Atrial Appendage'
GENE <- 'UST'
CELL <- 'Monocytes'

tis <- 'Thyroid'
GENE <- 'COMMD3' # 'DNAJC1'
CELL <- 'T cells follicular helper'

# df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)
df.infil <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

# infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
# df.pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt',data.table = F,stringsAsFactors = F)

EXPR.df <- df.geno[which(df.geno$Description==GENE),-c(1:2)]
EXPR.df <- as.data.frame(t(EXPR.df))
colnames(EXPR.df) <- 'GE'
EXPR.df$SAMP <- rownames(EXPR.df); rownames(EXPR.df) <- NULL
df.infil.sub <- subset(df.infil,SMTSD == tis)
dataf <- merge(df.infil.sub,EXPR.df,by.x='Input Sample',by.y='SAMP')

cor.test(dataf[,CELL],dataf$GE,use='p')







# phenoNum <- rownames(subset(infil_pheno,tissue==tis & cell==CELL))
# dataf <- df.pheno[,c("IID",paste0('pheno',phenoNum,'.2'))]



df.tmp <- data.frame(EXPR=as.numeric(df.geno[which(df.geno$GENE==GENE),-1]),CELL=df.infil[,CELL],ID=df.infil$ID)
cor.test(df.tmp[,1],df.tmp[,2])

g1 <- ggplot(df.tmp,aes(x=EXPR,y=CELL)) + geom_point(col='orange',pch=20) + geom_smooth(method='lm',se=F,col='black') +
  labs(x='STAM2 Expression',y='CD4+ T cells (CIBERSORT - Absolute)') + 
  theme_bw() + 
  theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank())

# save
scale <- 3
tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','Expr_x_Pheno','.png'),width=400*scale,height=400*scale,res=100*scale,units="px")
print(g1)
dev.off()


##############################

tis.nospace <- 'Artery_-_Tibial'
GENE <- 'KCTD10'
CELL <- 'MacrophageSum'
tis <- str_replace_all(tis.nospace,'_',' ')

f <- paste0('/Volumes/SeagateBackupPlusDrive/Elemento/GTEx_geno/subset_plink/gtex_all.filter.name.',tis.nospace,'.fam')
df.pheno <- fread(f,data.table = F,stringsAsFactors = F)
df.geno <- fread(paste0('/Volumes/SeagateBackupPlusDrive/Elemento/GTEx_GeneExpr_SplitByTissue/GTEx_v7_',tis.nospace,'_genexpr.txt'),data.table=F,stringsAsFactors = F)
df.infil <- fread('/Volumes/SeagateBackupPlusDrive/Elemento/CIBERSORT_out/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.infil <- subset(df.infil,SMTSD==tis)

df.tmp2 <- data.frame(EXPR=as.numeric(df.geno[which(df.geno$GENE==GENE),-1]),CELL=df.infil[,CELL],ID=df.infil$ID)
cor.test(df.tmp[,1],df.tmp[,2])

g2 <- ggplot(df.tmp2,aes(x=EXPR,y=CELL)) + geom_point(col='orange',pch=20) + geom_smooth(method='lm',se=F,col='black') +
  labs(x='KCTD10 Expression',y='Macrophages (CIBERSORT - Absolute)')  +
  theme_bw() + 
  theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank())

g2
library(cowplot)
# plot_grid(g1,g2)
scale <- 3
tiff(paste0('/Volumes/SeagateBackupPlusDrive/Elemento/120418/','Expr_x_Pheno','.png'),width=650*scale,height=500*scale,res=100*scale,units="px")
plot_grid(g2,g1,ncol=2)
dev.off()



##########################

library(stringr)
library(data.table)

GENE <- 'CUX1'
# GENE <- 'PTGS2'
# GENE <- 'IL1A'
# GENE <- 'MMP10'
# GENE <- 'F2RL1'
# GENE <- 'MMP9'
CELL <- 'Neutrophils'
tis.nospace <- 'Lung'
tis <- str_replace_all(tis.nospace,'_',' ')

# tis.nospace <- 'Skin_-_Sun_Exposed_(Lower_leg)'
# GENE <- 'STAM2'
# CELL <- 'CD4_Tcells'
# tis <- str_replace_all(tis.nospace,'_',' ')

tis.nospace <- 'Artery_-_Tibial'
GENE <- 'KCTD10'
CELL <- 'MacrophageSum'
tis <- str_replace_all(tis.nospace,'_',' ')

# df.geno <- fread('/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',data.table = F,stringsAsFactors = F)
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)
ind <- which(colnames(df.geno) %in% subset(df.attr,SMTSD==tis)$SAMPID)
df.geno.sub <- df.geno[,c(2,ind)]
df.infil <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.infil <- subset(df.infil,SMTSD==tis)
# which(colnames(df.geno.sub)[-1]!=df.infil$"Input Sample")
df.tmp <- data.frame(EXPR=as.numeric(df.geno.sub[which(df.geno.sub$Description==GENE),-1]),CELL=df.infil[,CELL],ID=df.infil$ID)
# cor.test(df.tmp[,1],df.tmp[,2])
summary(lm(df.tmp[,2]~df.tmp[,1]))$coef[2,4]
df.tmp$EXPR1=ceiling(df.tmp$EXPR / 30)
# df.tmp$EXPR2=ceiling(df.tmp$EXPR / 0.2)
fwrite(df.tmp,paste0('/athena/elementolab/scratch/anm2868/GTEx/',tis,'.',CELL,'.',GENE,'.txt'),col.names = T,row.names = F,sep='\t',quote=F)
df.tmp$EXPR1=ceiling(df.tmp$EXPR / 30)
table(df.tmp$EXPR1)
aggregate(df.tmp$CELL,list(df.tmp$EXPR1),median)


df.tmp2 <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/Lung.Neutrophils.CUX1.txt',stringsAsFactors = F,header = T,data.table = F)
# table(df.tmp2$EXPR1)
df.tmp2$EXPR1[df.tmp2$EXPR1>6] <- 6
summary(lm(df.tmp2$CELL~df.tmp2$EXPR1))
# g2 <- ggplot(df.tmp2,aes(x=EXPR1,y=CELL)) + geom_point(col='orange',pch=20) + geom_smooth(method='lm',se=F,col='black') +
#   labs(x='CUX1 Expression',y='Neutrophils (CIBERSORT - Absolute)')  +
#   theme_bw() + 
#   theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank())
g2 <- ggplot(df.tmp2,aes(x=as.factor(EXPR1),y=CELL)) +
  geom_boxplot(na.rm = T,outlier.shape = NA) +
  geom_jitter(col='steelblue',pch=20,width=0.1) + 
  labs(x=expression(paste(italic('CUX1'),' Expression (TPM)')),y='Neutrophil score (CIBERSORT - Absolute)')  +
  scale_x_discrete(label=c('10 - 15','15 - 20','25 - 30','30-40')) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank())

g2
library(cowplot)
# plot_grid(g1,g2)
scale <- 3
tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','Expr_x_Pheno','.png'),width=400*scale,height=400*scale,res=100*scale,units="px")
print(g2)
dev.off()

# tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','Expr_x_Pheno','.png'),width=650*scale,height=500*scale,res=100*scale,units="px")
# plot_grid(g2,g2,ncol=2)
# dev.off()


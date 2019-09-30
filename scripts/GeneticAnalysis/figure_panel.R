library(data.table)
library(ggplot2)
library(GenABEL)
library(cowplot)
library(gtools)

cluster_run <- FALSE

# 1: identify top 3 hits
# run on cluster:
if (cluster_run) {
  thres <- '5e-8'
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.ALL_EBM.CLUMPED_',thres,'.txt.clumped')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  x <- df[order(df$P),c('tissue','cell','SNP','P')][1:3,]
  infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
  i1 <- which(paste(infil_pheno$tissue,infil_pheno$cell)%in% paste(x$tissue[1],x$cell[1]))
  i2 <- which(paste(infil_pheno$tissue,infil_pheno$cell)%in% paste(x$tissue[2],x$cell[2]))
  i3 <- which(paste(infil_pheno$tissue,infil_pheno$cell)%in% paste('Esophagus - Muscularis','MastSum'))
  f <- '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt'
  df.pheno <- fread(f,data.table = F,stringsAsFactors = F)
  fwrite(df.pheno[,c('IID',paste0('pheno',c(i1,i2,i3),'.2'))],
         '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.top3.txt',
         col.names = T,row.names = F,quote=F,sep='\t',na='NA')
}

# initialize:
id <- list(); 
id[[1]] <- 166
id[[2]] <- 56
id[[3]] <- 81

# 2: plot qqplots for top 3 hits
df.qq <- list(); g.qq <- list()
# ggplot_title <- list('Lymphocytes in sigmoid colon samples','Monocytes in heart (atrial appendage) samples','Tfh cells in thyroid samples')
for (i in 1:3) {
  print(i)
  workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/'
  f <- paste0(workdir,'output/GeneticAnalysis/GWAS/','GTEx.pheno',id[[i]],'.ALL_EBM.txt')
  df <- fread(f,data.table = F,stringsAsFactors = F)
  pvector <- df$Pval_Brown
  n <- length(pvector)+1
  exp.x <- -log10((rank(pvector, ties.method="first")-.5)/n)
  pvalues <- -log10(pvector)
  
  df.qq[[i]] <- data.frame(exp=exp.x,obs=pvalues)
  
  # df.qq.full[[i]] <- data.frame(exp=exp.x,obs=pvalues)
  # df.qq[[i]] <- df.qq.full[[i]][sample(1:nrow(df.qq.full[[i]]),10000),]
  
  g.qq[[i]] <- ggplot(df.qq[[i]],aes(x=exp,y=obs)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
    theme_bw() + theme(panel.grid=element_blank(),plot.title = element_text(hjust=0.5)) + 
    scale_y_continuous(breaks=seq(0,max(pvalues),by=2)) +
    scale_x_continuous(breaks=seq(0,max(exp.x),by=2)) +
    labs(x = expression(Expected ~ ~-log[10](italic(p))), y = expression(Observed ~ ~-log[10](italic(p)))) #+
}

# 3: plot genotype phenotype plots for top 3 hits
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.top3.txt'
df.pheno <- fread(f,stringsAsFactors = F,data.table = F)

snp <- list(); rs <- list(); g.geno <- list(); phenoName <- list(); df.snp <- list(); pheno <- list(); df.geno <- list()
tis <- list('Thyroid','Colon - Sigmoid','Esophagus - Muscularis')
phenoName <- list("Tfh cells (trans resid)",
                  "Lymphocytes (trans resid)",
                  "Mast cells (trans resid)")
pheno <- as.list(paste0('pheno',id,'.2'))
snp <- list('10_22337752_A_G_b37','22_23961126_C_T_b37','17_77995143_A_G_b37')
rs <- list('rs6482199','rs56234965','rs9989443')
A1 <- list(); A2 <- list()
for (i in 1:3) {
  print(i)
  f <- paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/Genotype_Subset/gtex_all.filter.name.',snp[[i]],'.raw')
  df.snp[[i]] <- fread(f,data.table = F,stringsAsFactors = F)
  df.geno[[i]] <- merge(df.pheno[,c('IID',pheno[[i]])],df.snp[[i]][,c(2,grep(snp[[i]],colnames(df.snp[[i]])))],by='IID')
  df.geno[[i]] <- subset(df.geno[[i]],!(is.na(df.geno[[i]][,grep(snp[[i]],colnames(df.geno[[i]]))]) | (df.geno[[i]][,pheno[[i]]] == -9)))
  colnames(df.geno[[i]])[2] <- 'resid'
  colnames(df.geno[[i]])[ncol(df.geno[[i]])] <- 'genotype'
  A1[[i]] <- strsplit(snp[[i]],'_')[[1]][3]
  A2[[i]] <- strsplit(snp[[i]],'_')[[1]][4]
  g.geno[[i]] <- ggplot(df.geno[[i]],aes(x=as.factor(genotype),y=resid,fill=as.factor(genotype))) + geom_boxplot(na.rm = T,outlier.shape = NA) + geom_jitter(width=0.1,alpha=0.1) +
    labs(x=rs[[i]],y=phenoName[[i]]) + #,title=tis[[i]]) +
    theme_bw() +
    theme(legend.position="none",plot.title = element_text(hjust=0.5)) +
    theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank()) +
    scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
    scale_x_discrete(labels=c(paste0(A1[[i]],A1[[i]]),paste0(A1[[i]],A2[[i]]),paste0(A2[[i]],A2[[i]])))
  
}
library(cowplot)
# plot_grid(g.geno[[1]],g.geno[[2]],g.geno[[3]],ncol=3)

####################################################################

# 4: plot expression by phenotype plots

# for running on cluster:
if (cluster_run) {
  df.exp <- fread('/athena/elementolab/scratch/anm2868/GTEx/EXPRESSION/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',data.table = F,stringsAsFactors = F)
  GENE <- list('COMMD3','C22orf43','CCDC40');
  EXPR.df <- df.exp[which(df.exp$Description %in% GENE),-1]
  rownames(EXPR.df) <- EXPR.df$Description; EXPR.df <- EXPR.df[,-1]
  EXPR.df <- as.data.frame(t(EXPR.df))
  EXPR.df$SAMP <- rownames(EXPR.df); rownames(EXPR.df) <- NULL
  paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
  EXPR.df$SAMP2 <- sapply(strsplit(EXPR.df$SAMP,"-"),paste.s)
  fwrite(EXPR.df,
         '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/expression_data.txt',
         col.names = T,row.names = F,quote=F,sep='\t',na='NA')
}

#
EXPR.df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/expression_data.txt',data.table = F,stringsAsFactors = F)
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.top3.txt'
df.pheno <- fread(f,stringsAsFactors = F,data.table = F)
df.infil <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

CELL <- list('T cells follicular helper','Lymph_Sum','MastSum')
tis <- list('Thyroid','Colon - Sigmoid','Esophagus - Muscularis')
GENE <- list('COMMD3','C22orf43','CCDC40');
phenoName <- list("Tfh cells (transformed)",
                  "Lymphocytes (transformed)",
                  "Mast cells (transformed)")
pheno <- as.list(paste0('pheno',id,'.2'))

g.expr1 <- list(); g.expr2 <- list(); res <- list()
for (i in 1:3) {
  print(i)
  
  # df.pheno.sub <- subset(df.pheno,df.pheno[,pheno[[i]]]!=-9)[,c('IID',pheno[[i]])]
  
  # if want to use original infiltration data:
  X <- subset(df.infil,SMTSD==tis[[i]])
  df.pheno.sub <- X[,c('ID',CELL[[i]])]; colnames(df.pheno.sub) <- c('IID',pheno[[i]])
  
  full_id <- X[,'Input Sample']
  EXPR.df.sub <- subset(EXPR.df,SAMP %in% full_id)
  
  df.mg <- merge(EXPR.df.sub[,c('SAMP2',GENE[[i]])],df.pheno.sub[,c('IID',pheno[[i]])],by.x='SAMP2',by.y='IID')
  
  # plot:
  df.mg$EXPR_BIN <- quantcut(df.mg[,GENE[[i]]],3)
  df.mg$PHENOTYPE <- df.mg[,pheno[[i]]]
  geneName <- GENE[[i]]
  x_axis_label <- bquote(italic(.(geneName))~'Expression (TPM)')
  g.expr2[[i]] <- ggplot(df.mg,aes(x=as.factor(EXPR_BIN),
                                  # y=(PHENOTYPE))) +
                                  y=rntransform(PHENOTYPE))) +
    geom_boxplot(na.rm = T,outlier.shape = NA,width=0.5) +
    geom_jitter(col='steelblue',pch=20,width=0.1,alpha=0.25) +
    labs(x=x_axis_label,y=phenoName[[i]])  +
    theme_bw() + 
    theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank())
  
  res[[i]] <- cor.test(df.mg[,GENE[[i]]],rntransform(df.mg$PHENOTYPE))
  
}

figure_panel <- plot_grid(g.qq[[1]],g.geno[[1]],g.expr2[[1]],
                           g.qq[[2]],g.geno[[2]],g.expr2[[2]],
                           g.qq[[3]],g.geno[[3]],g.expr2[[3]],
                           ncol=3,rel_widths = c(0.85,1,1))
# figure_panel
tiff('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/figure_panel.png',width=4500,height=5000,res=400,units="px")
print(figure_panel)
dev.off()

print(res)

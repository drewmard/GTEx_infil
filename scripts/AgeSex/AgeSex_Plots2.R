# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: plot of age & sex influences over infiltration patterns

# load & initialize
library(ggplot2)
library(GenABEL)
library(data.table)
workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/'
df.abs <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
df.rel <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
df.xcell <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt'),data.table = F,stringsAsFactors = F)

# lymphocytes in bresat & thyroid tissue
tissue <- list('Breast - Mammary Tissue','Thyroid')
g.plots <- list(); df.mg.melt <- list()
nplot <- 2
for (i in 1:nplot) {
  tis <- tissue[[i]]
  df.mg <- subset(df.abs,SMTSD==tis)[,c('ID','SEX','Lymph_Sum')]
  colnames(df.mg)[3] <- c('CIB-Abs')
  # for (j in 3:3) df.mg[,j] <- rntransform(df.mg[,j])
  df.mg.melt[[i]] <- (melt(df.mg[,-1],id='SEX'))
  g.plots[[i]] <- ggplot(df.mg.melt[[i]],aes(x=as.factor(SEX),y=value,fill=as.factor(SEX))) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.1) +
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) + theme(legend.position='none',panel.grid = element_blank()) +
    # labs(x='Sex',y='Lymphocytes (transformed scores)',title = tis) + scale_x_discrete(labels=c('Male','Female')); #g.plots[[i]]
    labs(x='Sex',y='Lymphocytes score',title = tis) + scale_x_discrete(labels=c('Male','Female')); #g.plots[[i]]
  if (i!=nplot) {g.plots[[i]] <- g.plots[[i]] + theme(legend.position='none')}
}
library(cowplot)
plot_grid(g.plots[[1]],g.plots[[2]],ncol=2,rel_widths = c(0.5,0.5))
scale <- 3
tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/BreastThyroid.png",width=800*scale,height=600*scale,res=100*scale,units="px")
plot_grid(g.plots[[1]],g.plots[[2]],ncol=2,rel_widths = c(0.5,0.5))
dev.off()


# tibial-age associations 

cellTypes.df <- data.frame( 
  ciber=c('T cells CD8','T cells CD4 naive','CD4_memory','Neutrophils','MacrophageSum',
          'Bcellsum','NK_Sum','DendriticSum','MastSum','Myeloid_Sum',
          'T cells follicular helper','T cells regulatory (Tregs)','T cells gamma delta',
          'Monocytes','Eosinophils','Lymph_Sum'),
  xcell=c('CD8Sum','CD4+ naive T-cells','CD4_memory','Neutrophils','MacrophageSum',
          'Bcellsum','NK cells','DendriticSum','Mast cells','Myeloid_Sum',
          'Th_Sum','Tregs','Tgd cells',
          'Monocytes','Eosinophils','Lymph_Sum'),
  stringsAsFactors = F)

tissue <- list('Artery - Tibial','Nerve - Tibial')
g.plots <- list(); df.mg.melt <- list()
nplot <- 2
celltypes <- c('T cells CD8','CD4_memory','MacrophageSum','MastSum','Monocytes','Myeloid_Sum','Lymph_Sum')
celltypes.col <- c('CD8+ T cells','CD4+ memory T cells','Macrophages','Mast cells','Monocytes','Myeloid cells','Lymphocytes')
for (i in 1:nplot) {
  tis <- tissue[[i]]
  if (i==1) {celltypes2 <- celltypes[-1];celltypes2.col <- celltypes.col[-1]} else {celltypes2 <- celltypes; celltypes2.col <- celltypes.col}
  df.mg <- subset(df.abs,SMTSD %in% tis)[,c('ID','AGE','SMTSD',celltypes2)]
  df.mg.melt[[i]] <- (melt(df.mg[,-1],id=c('AGE','SMTSD')))
  g.plots[[i]] <-
    ggplot(df.mg.melt[[i]],aes(x=as.factor(variable),y=value,fill=as.factor(AGE))) + 
    geom_boxplot(outlier.alpha=0.1) +
    scale_x_discrete(label=celltypes2.col) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) +#+ theme(legend.position='none',panel.grid = element_blank()) +
    labs(x='Cell type',y='Scores',fill='Age',title=tis) #+ scale_x_discrete(labels=c('Male','Female')); #g.plots[[i]]
  #if (i!=nplot) {g.plots[[i]] <- g.plots[[i]] + theme(legend.position='none')}
}
library(cowplot)
plot_grid(g.plots[[1]],g.plots[[2]],ncol=1)#,rel_widths = c(0.5,0.6))
tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/Tibial.png",width=800*scale,height=400*scale,res=60*scale,units="px")
plot_grid(g.plots[[1]],g.plots[[2]],ncol=1,rel_widths = c(0.5,0.5))
dev.off()






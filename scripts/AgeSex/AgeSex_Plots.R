# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: plot of age & sex influences over infiltration patterns

library(ggplot2)
workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/'
df.abs <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
df.rel <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
df.xcell <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt'),data.table = F,stringsAsFactors = F)

tissue <- list('Breast - Mammary Tissue','Thyroid')
g.plots <- list(); df.mg.melt <- list()
nplot <- 2
for (i in 1:nplot) {
  tis <- tissue[[i]]
  df.mg <- merge(
    merge(
      subset(df.abs,SMTSD==tis)[,c('ID','SEX','Lymph_Sum')],
      subset(df.rel,SMTSD==tis)[,c('ID','Lymph_Sum')],by='ID'),
    subset(df.xcell,SMTSD==tis)[,c('IID','Lymph_Sum')],
    by.x='ID',by.y='IID')
  colnames(df.mg)[3:5] <- c('CIB-Abs','CIB-Rel','xCell')
  for (j in 3:5) df.mg[,j] <- rntransform(df.mg[,j])
  df.mg.melt[[i]] <- (melt(df.mg[,-1],id='SEX'))
  g.plots[[i]] <- ggplot(df.mg.melt[[i]],aes(x=as.factor(SEX),y=value,fill=as.factor(variable))) + geom_boxplot(outlier.shape=NA) +
    #geom_jitter(width=0.1) + 
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) +#+ theme(legend.position='none',panel.grid = element_blank()) +
    labs(x='Sex',y='Lymphocytes (transformed scores)',title = tis,fill='Method') + scale_x_discrete(labels=c('Male','Female')); #g.plots[[i]]
  if (i!=nplot) {g.plots[[i]] <- g.plots[[i]] + theme(legend.position='none')}
}
library(cowplot)
plot_grid(g.plots[[1]],g.plots[[2]],ncol=2,rel_widths = c(0.5,0.6))

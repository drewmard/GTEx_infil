library(data.table)
library(ggplot2)

df.ciber.abs <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table=F)
df.ciber.rel <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table=F,stringsAsFactors = F)
df.xcell <- fread('~/Documents/Research/GTEx/Infiltration/xCell_TPM_TBT_mod.txt',data.table=F,stringsAsFactors = F)

# initialize:
dataf <- data.frame(ID=unique(df.ciber.abs$ID))
s <- data.frame(table(df.ciber.rel$SMTSD))

df.ciber.abs <- subset(df.ciber.abs,SMTSD %in% s[,1])
colnames(df.ciber.abs)[colnames(df.ciber.abs)=='CD4_Tcells'] <- 'CD4Sum'
colnames(df.ciber.abs)[colnames(df.ciber.abs)=='T cells CD8'] <- 'CD8Sum'

df.ciber.rel <- subset(df.ciber.rel,SMTSD %in% s[,1])
colnames(df.ciber.rel)[colnames(df.ciber.rel)=='CD4_Tcells'] <- 'CD4Sum'
colnames(df.ciber.rel)[colnames(df.ciber.rel)=='T cells CD8'] <- 'CD8Sum'

df.xcell <-  subset(df.xcell,SMTSD %in% s[,1])

for (i in 1:4) {
  
  if (i==1) {
    cell <- 'MacrophageSum'; cell_name <- 'Macrophages'
  } else if (i==2) {
    cell <- 'Neutrophils'; cell_name <- cell
  } else if (i==3) {
    cell <- 'CD4Sum'; cell_name <- 'CD4+ T cells'
  } else if (i==4) {
    cell <- 'CD8Sum'; cell_name <- 'CD8+ T cells'
  }

  # library(stringr)
  # s[,1] <- str_replace_all(s[,1],'_',' ') 

  df.ciber.rel2 <- merge(df.ciber.rel,df.ciber.abs[,c('Input Sample',cell)],by = 'Input Sample')
  df.xcell2 <- merge(df.xcell,df.ciber.abs[,c('Input Sample',cell)],by.x = 'SAMP',by.y='Input Sample')
  
  g1 <- ggplot(df.ciber.rel2,aes(x=reorder(SMTSD,df.ciber.rel2[,paste0(cell,'.y')],FUN=median),y=df.ciber.rel2[,paste0(cell,'.x')])) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
    labs(x='Sample Type',y=cell_name);
  
  g2 <- ggplot(df.ciber.abs,aes(x=reorder(SMTSD,df.ciber.abs[,cell],FUN=median),y=df.ciber.abs[,cell])) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
    labs(x='Sample Type',y=cell_name);
  
  g3 <- ggplot(df.xcell2,aes(x=reorder(SMTSD,df.xcell2[,paste0(cell,'.y')],FUN=median),y=(df.xcell2[,paste0(cell,'.x')])^1)) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 80, hjust = 1,size = rel(1)),panel.grid = element_blank()) +
    labs(x='Sample Type',y=cell_name);
  
  # to plot together:
  # library(cowplot)
  # plot_grid(g1,g2,g3,ncol=1)
  
  # to save plots:
  x <- 3
  tiff(paste0('~/Documents/Research/GTEx/Infiltration/',cell,'_rel.png'),width=1000*x,height=400*x,res=100*x,units="px")
  print(g1)
  dev.off()
  tiff(paste0('~/Documents/Research/GTEx/Infiltration/',cell,'_abs.png'),width=1000*x,height=400*x,res=100*x,units="px")
  print(g2)
  dev.off()
  tiff(paste0('~/Documents/Research/GTEx/Infiltration/',cell,'_xcell.png'),width=1000*x,height=400*x,res=100*x,units="px")
  print(g3)
  dev.off()

}

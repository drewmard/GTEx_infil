# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: generate boxplots to compare cell type scores across tissues and between deconvolution methods

# load
library(data.table)
library(ggplot2)
library('cowplot')
# cell type of interest:
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
celltype.label.vec <- c('CD8+ T cells','CD4+ naive T cells','CD4+ memory T cells','Neutrophils','Macrophages',
                        'B cells','NK cells','Dendritic cells','Mast cells','Myeloid cells',
                        'Helper T cells','Tregs','Gamma delta T cells',
                        'Monocytes','Eosinophils','Lymphocytes')

# load infiltration data:
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
x <- as.data.frame(table(df$SMTSD)); tis <- subset(x,Freq>=70)$Var1; tis <- tis[!(tis %in% c('Cells - Transformed fibroblasts','Cells - EBV-transformed lymphocytes'))]
df <- subset(df, SMTSD %in% tis)
ciber.rel <- df[,c('Input Sample','SMTSD',cellTypes.df$ciber)]
colnames(ciber.rel)[-c(1:2)] <- paste0(cellTypes.df$ciber,'.CIB_Rel')

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
ciber.abs <- df[,c('Input Sample',cellTypes.df$ciber)]
colnames(ciber.abs)[-1] <- paste0(cellTypes.df$ciber,'.CIB_Abs')

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
xcell <- df[,c('SAMP',cellTypes.df$xcell)]
colnames(xcell)[-1] <- paste0(cellTypes.df$ciber,'.xCell')

df <- merge(merge(ciber.rel,ciber.abs,by='Input Sample'),xcell,by.x='Input Sample',by.y='SAMP')

for (i in 1:nrow(cellTypes.df)) {
  
  celltype <- cellTypes.df$ciber[i]
  celltype.xcell <- subset(cellTypes.df,ciber==celltype)$xcell
  celltype.label <- celltype.label.vec[i]
  print(paste0('cell type: ',celltype.label))
  
  # create plots
  g1 <- ggplot(df,aes(x=reorder(SMTSD,df[,paste0(celltype,'.CIB_Abs')],FUN=median),y=df[,paste0(celltype,'.CIB_Rel')])) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x='Tissue Type',y=celltype.label);# g1
  
  g2 <- ggplot(df,aes(x=reorder(SMTSD,df[,paste0(celltype,'.CIB_Abs')],FUN=median),y=df[,paste0(celltype,'.CIB_Abs')])) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x='Tissue Type',y=celltype.label);# g2
  
  g3 <- ggplot(df,aes(x=reorder(SMTSD,df[,paste0(celltype,'.CIB_Abs')],FUN=median),y=df[,paste0(celltype,'.xCell')])) + 
    geom_point(col='orange',alpha=0.6) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
    labs(x='Tissue Type',y=celltype.label);# g3
  
  x <- 3
  dir.create('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/pantissue_plots',showWarnings = F)
  prefix <- paste0('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/pantissue_plots/',celltype.label,'_boxplot')
  
  # single plots:
  # tiff(paste0(prefix,'_rel.png'),width=1000*x,height=400*x,res=100*x,units="px")
  # print(g1)
  # dev.off()
  # tiff(paste0(prefix,'_abs.png'),width=1000*x,height=400*x,res=100*x,units="px")
  # print(g2)
  # dev.off()
  # tiff(paste0(prefix,'_xcell.png'),width=1000*x,height=400*x,res=100*x,units="px")
  # print(g3)
  # dev.off()
  
  # combined plots:
  # tiff(paste0(prefix,'_All3.png'),width=750*x,height=1050*x,res=100*x,units="px")
  tiff(paste0(prefix,'_All3.png'),width=750*x,height=900*x,res=100*x,units="px")
  print(plot_grid(g1,g2,g3,ncol=1,rel_heights = c(0.2,0.2,0.35)))
  dev.off()
  
}


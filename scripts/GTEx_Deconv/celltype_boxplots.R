library(data.table)
library(ggplot2)

# choose cell type from cellTypes.df$ciber:
celltype <- "MacrophageSum"
# choose y axis label
celltype.label <- 'Macrophages'

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
celltype.xcell <- subset(cellTypes.df,ciber==celltype)$xcell

# load infiltration data:
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
x <- as.data.frame(table(df$SMTSD)); tis <- subset(x,Freq>=70)$Var1; tis <- tis[!(tis %in% c('Cells - Transformed fibroblasts','Cells - EBV-transformed lymphocytes'))]
df <- subset(df, SMTSD %in% tis)
ciber.rel <- df[,c('SMTSD',cellTypes.df$ciber)]
ciber.rel <- melt(ciber.rel,id='SMTSD')
ciber.rel$method<-'Cib-Rel'

# save:
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
ciber.abs <- df[,c('SMTSD',cellTypes.df$ciber)]
ciber.abs <- melt(ciber.abs,id='SMTSD')
ciber.abs$method<-'Cib-Abs'

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
xcell <- df[,c('SMTSD',cellTypes.df$xcell)]
xcell <- melt(xcell,id='SMTSD')
xcell$method<-'xCell'

# merge all 3 methods
df <- rbind(rbind(ciber.rel,ciber.abs),xcell)
df$method <- factor(df$method,levels = c('Cib-Rel','Cib-Abs','xCell'))

df.sub <- subset(df,variable==celltype)
df.sub2 <- subset(df.sub,method=='Cib-Abs')

g1 <- ggplot(subset(df,method=='Cib-Rel' & variable==celltype),aes(x=reorder(SMTSD,subset(df,method=='Cib-Abs' & variable==celltype)$value,FUN=median),y=value)) + 
  geom_point(col='orange',alpha=0.6) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
  labs(x='Tissue Type',y=celltype.label);# g1

g2 <- ggplot(subset(df,method=='Cib-Abs' & variable==celltype),aes(x=reorder(SMTSD,subset(df,method=='Cib-Abs' & variable==celltype)$value,FUN=median),y=value)) + 
  geom_point(col='orange',alpha=0.6) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
  labs(x='Tissue Type',y=celltype.label);# g2

g3 <- ggplot(subset(df,method=='xCell' & variable==celltype),aes(x=reorder(SMTSD,subset(df,method=='Cib-Abs' & variable==celltype)$value,FUN=median),y=value)) + 
  geom_point(col='orange',alpha=0.6) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),panel.grid = element_blank()) +
  labs(x='Tissue Type',y=celltype.label);# g3

x <- 3
prefix <- paste0('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/',celltype.label,'_boxplot')
tiff(paste0(prefix,'_rel.png'),width=1000*x,height=400*x,res=100*x,units="px")
print(g1)
dev.off()
tiff(paste0(prefix,'_abs.png'),width=1000*x,height=400*x,res=100*x,units="px")
print(g2)
dev.off()
tiff(paste0(prefix,'_xcell.png'),width=1000*x,height=400*x,res=100*x,units="px")
print(g3)
dev.off()



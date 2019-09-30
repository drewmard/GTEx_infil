# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script for building clustered heatmaps showing tissue variability across different cell types

library(data.table)
library(pheatmap)
library(matrixStats)

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
cellTypes.df <- subset(cellTypes.df,!(ciber %in% c('Myeloid_Sum','Lymph_Sum')))

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
x <- as.data.frame(table(df$SMTSD)); tis <- subset(x,Freq>=70)$Var1; tis <- tis[!(tis %in% c('Cells - Transformed fibroblasts','Cells - EBV-transformed lymphocytes'))]
df <- subset(df, SMTSD %in% tis)
ciber.rel <- df[,c('Input Sample','SMTSD',cellTypes.df$ciber)]

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
ciber.abs <- df[,c('Input Sample',cellTypes.df$ciber)]


df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
df <- subset(df, SMTSD %in% tis)
xcell <- df[,c('SAMP',cellTypes.df$xcell)]

# merge all 3 methods
df <- merge(merge(ciber.rel,ciber.abs,by='Input Sample'),xcell,by.x='Input Sample',by.y='SAMP')

tissueMdns <- matrix(, nrow = length(tis), ncol = nrow(cellTypes.df)*3)
for (i in 1:length(tis)) {
  x <- subset(df[,-1],SMTSD %in% tis[i])[,-1]
  tissueMdns[i,] <- colMedians(as.matrix(x))
  
}
dataf <- as.numeric(tissueMdns[,1:nrow(cellTypes.df)])
dataf[which(dataf > 0.3)] <- 0.3
y<-matrix(scale(dataf),nrow=length(tis),ncol=nrow(cellTypes.df))
tissueMdns[,1:nrow(cellTypes.df)] <- y

dataf <- as.numeric(tissueMdns[,(nrow(cellTypes.df)+1):(2*nrow(cellTypes.df))])
dataf[which(dataf > 0.6)] <- 0.6
y<-matrix(scale(dataf),nrow=length(tis),ncol=nrow(cellTypes.df))
tissueMdns[,(nrow(cellTypes.df)+1):(2*nrow(cellTypes.df))] <- y

dataf <- as.numeric(tissueMdns[,(2*nrow(cellTypes.df)+1):(3*nrow(cellTypes.df))])
# dataf[which(dataf > 0.1)] <- 0.1
dataf[which(dataf > 0.05)] <- 0.05
y<-matrix(scale(dataf),nrow=length(tis),ncol=nrow(cellTypes.df))
tissueMdns[,(2*nrow(cellTypes.df)+1):(3*nrow(cellTypes.df))] <- y

tissueMdns <- t(tissueMdns)
colnames(tissueMdns) <- tis
rownames(tissueMdns) <- paste0(rep(c('CD8 T cells','CD4 naive T cells','CD4 memory T cells',
        'Neutrophils','Macrophages','B cells',
        'NK cells','Dendritic cells','Mast cells',
        'Helper T cells','Tregs','Tgd cells',
        'Monocytes','Eosinophils'),3),
       rep(c(' (Cib-Rel)',' (Cib-Abs)',' (xCell)'),each=14)
)

i <- rep(seq(1,nrow(cellTypes.df)),each=3)
i[seq(2,nrow(cellTypes.df)*3,by=3)] <- i[seq(2,nrow(cellTypes.df)*3,by=3)] + nrow(cellTypes.df)
i[seq(3,nrow(cellTypes.df)*3,by=3)] <- i[seq(3,nrow(cellTypes.df)*3,by=3)] + nrow(cellTypes.df)*2
tissueMdns.ordered <- tissueMdns[i,]
# p <- pheatmap(tissueMdns.ordered, cluster_rows = FALSE,clustering_distance_cols = "euclidean", 
#               clustering_method = "complete",cex = 0.8, angle_col = 45, 
#               xlab = "Tissues",#breaks = myBreaks,
#               cellheight=5,cellwidth=10,
#               annotation_colors = rep('black','orange','yellow'),
#               legend=F,
#               fontsize = 8); p

p.rel <- pheatmap(tissueMdns[1:nrow(cellTypes.df),], cluster_rows = FALSE,clustering_distance_cols = "euclidean", 
              clustering_method = "complete",cex = 0.8, angle_col = 45, 
              xlab = "Tissues",#breaks = myBreaks,
              cellheight=5,cellwidth=10,
              annotation_colors = rep('black','orange','yellow'),
              legend=F,
              fontsize = 8);# p.rel
p.abs <- pheatmap(tissueMdns[(nrow(cellTypes.df)+1):(nrow(cellTypes.df)*2),], cluster_rows = FALSE,clustering_distance_cols = "euclidean", 
                  clustering_method = "complete",cex = 0.8, angle_col = 45, 
                  xlab = "Tissues",#breaks = myBreaks,
                  cellheight=5,cellwidth=10,
                  annotation_colors = rep('black','orange','yellow'),
                  legend=T,
                  fontsize = 8); p.abs
p.xcell <- pheatmap(tissueMdns[(2*nrow(cellTypes.df)+1):(nrow(cellTypes.df)*3),], cluster_rows = FALSE,clustering_distance_cols = "euclidean", 
                  clustering_method = "complete",cex = 0.8, angle_col = 45, 
                  xlab = "Tissues",#breaks = myBreaks,
                  cellheight=5,cellwidth=10,
                  annotation_colors = rep('black','orange','yellow'),
                  legend=F,
                  fontsize = 8); #p.xcell

# scale <- 3
# tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/HeatMap.png",width=900*scale,height=600*scale,res=100*scale,units="px")
# print(p)
# dev.off()

scale <- 3
tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/HeatMap.xCell.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(p.xcell)
dev.off()

scale <- 3
tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/HeatMap.Cib-abs.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(p.abs)
dev.off()

scale <- 3
tiff("~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/HeatMap.Cib-rel.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(p.rel)
dev.off()

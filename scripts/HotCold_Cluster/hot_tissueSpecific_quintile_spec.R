library(ggplot2)
library(data.table)
library(pheatmap)
library(matrixStats)
library(cowplot)
k=7; i=k
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_phenotypes.txt',data.table = F,stringsAsFactors = F)
df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/FullHotColdAssignments_Correct.txt',data.table = F,stringsAsFactors = F)
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ConsensusQuintileAssignments-FullSecond.txt',data.table = F,stringsAsFactors = F); colnames(df)[which(colnames(df)=='ID')] <- 'IID'
df$IID <- as.character(lapply(strsplit(df$Sample,"\\."),function(x) paste(x[1:2],collapse = "-")))
# df$Top <- as.numeric(df$Quintile%in%c('Top','Second'))
# df$Top <- as.numeric(df$Quintile%in%c('Top'))
df$Top <- as.numeric(df$Infiltration %in% c('Hot'))
CellType.unique <- unique(df$CellType)
CellType.names <- c("CD4 memory T cells",'Macrophages','Mast cells','Myeloid cells','Monocytes','Lymphocytes','CD8 T cells',
                    'B cells','Dendritic cells','Helper T cells','Neutrophils','NK cells')
df.sub <- subset(df,CellType==CellType.unique[k])
x <- aggregate(df.sub$IID,list(df.sub$IID),length)
# df.sub <- subset(df.sub,x >= 8)
df.sub <- df.sub[,c('IID','Tissue','Top')]
library(reshape2)
df.sub.2 <- dcast(df.sub,IID~Tissue)
dataf <- df.sub.2[,c('Whole Blood','Lung')]
dataf <- dataf[apply(dataf,1,function(x) sum(is.na(x))==0),]
f <- paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/TissueHotHeatMap.',CellType.unique[i],'WB.Lung.consensus_clustering.png')
x=3;pheatmap(t(dataf),cluster_cols = T,cluster_rows=T,show_colnames = F,
         legend=F,treeheight_row = 50, treeheight_col = 0,border_color='black',
         cellwidth = rel(1.2*x),cellheight=rel(4*x),filename = f)


table(dataf)

# library(grid)
# p <- pheatmap(df.sub.5[,1:50],cluster_cols = T,cluster_rows=T,show_colnames = F,
#               legend=F,treeheight_row = 50, treeheight_col = 0,border_color='black',
#               cellwidth = rel(1.2*x),cellheight=rel(4*x),filename=f)
# my_gtable = p$gtable
# my_gtable$grobs[[3]]$gp=gpar(col="grey30",fontfamily='Helvetica',fontsize=7)
# p$gtable <- my_gtable
# print(p)



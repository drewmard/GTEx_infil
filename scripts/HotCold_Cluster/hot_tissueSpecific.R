# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: assess "hot" clusters across tissues within individuals


library(ggplot2)
library(data.table)
library(pheatmap)
library(matrixStats)
library(cowplot)

df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_phenotypes.txt',data.table = F,stringsAsFactors = F)
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ExpandedHotCold.txt',data.table = F,stringsAsFactors = F); colnames(df)[which(colnames(df)=='ID')] <- 'IID'
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ConsensusHotColdAssignments.txt',data.table = F,stringsAsFactors = F)
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/FullHotColdAssignments.txt',data.table = F,stringsAsFactors = F)
df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/FullHotColdAssignments_Correct.txt',data.table = F,stringsAsFactors = F)

df$IID <- as.character(lapply(strsplit(df$Sample,"\\."),function(x) paste(x[1:2],collapse = "-")))
df$Top <- as.numeric(df$Infiltration=='Hot')
df$Bottom <- as.numeric(df$Infiltration=='Cold')
df$Intermediate <- as.numeric(df$Infiltration=='Intermediate')
CellType.unique <- unique(df$CellType)
CellType.names <- c("CD4 memory T cells",'Macrophages','Mast cells','Myeloid cells','Monocytes','Lymphocytes','CD8 T cells',
                    'B cells','Dendritic cells','Helper T cells','Neutrophils','NK cells')

# 1 #############
# create density plots (supplement)
g <- list();dataf <- list();val <- list();val2 <- list()
for (i in 1:7) {
  df.sub <- subset(df,CellType==CellType.unique[i])
  tab <- aggregate(df.sub$Top,by=list(df.sub$IID),mean)
  colnames(tab) <- c('IID','Mean')
  tab2 <- aggregate(df.sub$Top,by=list(df.sub$IID),length)
  colnames(tab2) <- c('IID','N')
  tab <- merge(tab,tab2,by='IID')
  tab <- subset(tab,Mean>0 & N >= 8)
  dataf[[i]] <- as.data.frame(tab)
  val[[i]] <- median(dataf[[i]]$Mean)
  val2[[i]] <- median(dataf[[i]]$N)
  g[[i]] <- ggplot(dataf[[i]],aes(x=Mean)) + 
    geom_density(fill='orange3') + 
    geom_vline(xintercept = val[[i]],col='red',lty='dashed') +
    labs(x='Proportion of "hot" tissue samples',y='Number of individuals',title=paste0(CellType.names[i],' (N=',length(unique(df.sub$Tissue)),'; n=',val2[[i]],')'))+ 
    theme_bw() + scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),limits = c(-0.1,1.1)) +
    theme(plot.title=element_text(hjust=0.5))
}
library(cowplot)
plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],
          g[[5]],g[[6]],g[[7]],nrow=2)
mean(as.numeric(val))
as.numeric(lapply(dataf,function(x) {median(x$N)}))


# 2 #############
# create bar plots (supplement)
g <- list();dataf <- list();val <- list();val2 <- list()
for (i in 1:7) {
  df.sub <- subset(df,CellType==CellType.unique[i]); #unique(df.sub$Tissue)
  tab <- aggregate(df.sub$Top,by=list(df.sub$IID),sum)
  tab2 <- aggregate(df.sub$Top,by=list(df.sub$IID),length)
  colnames(tab) <- c('IID','Sum')
  colnames(tab2) <- c('IID','N')
  tab <- merge(tab,tab2,by='IID')
  tab <- subset(tab,Sum>0 & N >= 8)
  dataf[[i]] <- as.data.frame(tab)
  val[[i]] <- median(dataf[[i]]$Sum)
  val2[[i]] <- median(dataf[[i]]$N)
  g[[i]] <- ggplot(dataf[[i]],aes(x=as.integer(Sum))) + 
    geom_bar(fill='orange3') +
    labs(x='Number of "hot" tissue samples',y='Number of individuals',title=paste0(CellType.unique[i],' (N=',length(unique(df.sub$Tissue)),'; n=',val2[[i]],')'))+ 
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5)) +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
}

plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],
          g[[5]],g[[6]],g[[7]],nrow=2)
mean(as.numeric(val))
as.numeric(lapply(dataf,function(x) {median(x$N)}))
as.numeric(val)

# 3 #############

# create clustered heatmaps
for (i in 1:7) {
  df.sub <- subset(df,CellType==CellType.unique[i])
  x <- aggregate(df.sub$IID,list(df.sub$IID),length)
  df.sub <- subset(df.sub,x >= 8)
  df.sub <- df.sub[,c('IID','Tissue','Top')]
  library(reshape2)
  df.sub.2 <- dcast(df.sub,IID~Tissue)
  df.sub.3 <- melt(df.sub.2,id.vars = 'IID',na.rm = F)
  colnames(df.sub.3) <- c('IID','Tissue','Top')
  df.sub.3$Top[which(is.na(df.sub.3$Top))] <- 0
  df.sub.4 <- dcast(df.sub.3,IID~Tissue)
  rownames(df.sub.4) <- df.sub.4[,1]; df.sub.4 <- df.sub.4[,-1]
  df.sub.5 <- as.data.frame(t(df.sub.4))
  if (i==2) {fwrite(df.sub.5,'/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/Macrophage.3b.sourcedata.txt',quote = F,na='NA',sep = '\t',row.names = T,col.names = T)}
  x <- apply(df.sub.5,2,sum); df.sub.5 <- df.sub.5[,-which(x==0)]
  f <- paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/TissueHotHeatMap.',CellType.unique[i],'.png')
  x <-3; p <- pheatmap(df.sub.5,cluster_cols = T,cluster_rows=T,show_colnames = F,
                       legend=F,treeheight_row = 50, treeheight_col = 0,border_color='black',
                       cellwidth = rel(1.2*x),cellheight=rel(4*x),filename=f)
  library(grid)
  my_gtable = p$gtable
  my_gtable$grobs[[3]]$gp=gpar(col="grey30",fontfamily='Helvetica',fontsize=7)
  p$gtable <- my_gtable
  print(p)
  
}

######################################################
######################################################
######################################################

library(ggplot2)
library(data.table)
library(pheatmap)
library(matrixStats)
library(cowplot)
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_phenotypes.txt',data.table = F,stringsAsFactors = F)

df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/FullHotColdAssignments_Correct.txt',data.table = F,stringsAsFactors = F)
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ExpandedHotCold.txt',data.table = F,stringsAsFactors = F); colnames(df)[which(colnames(df)=='ID')] <- 'IID'
# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ConsensusHotColdAssignments.txt',data.table = F,stringsAsFactors = F)

df$IID <- as.character(lapply(strsplit(df$Sample,"\\."),function(x) paste(x[1:2],collapse = "-")))
df$Top <- as.numeric(df$Infiltration=='Hot')
df$Bottom <- as.numeric(df$Infiltration=='Cold')
df$Intermediate <- as.numeric(df$Infiltration=='Intermediate')
CellType.unique <- unique(df$CellType)
CellType.names <- c("CD4 memory T cells",'Macrophages','Mast cells','Myeloid cells','Monocytes','Lymphocytes','CD8 T cells',
                    'B cells','Dendritic cells','Helper T cells','Neutrophils','NK cells')


df.save.list <- list()
# for (k in 1:7) {
for (k in 1:7) {
  df.sub <- subset(df,CellType==CellType.unique[k])
  # x <- aggregate(df.sub$IID,list(df.sub$IID),length)
  # df.sub <- subset(df.sub,x >= 8)
  df.sub <- df.sub[,c('IID','Tissue','Top')]
  library(reshape2)
  df.sub.2 <- dcast(df.sub,IID~Tissue)
  df.sub.3 <- melt(df.sub.2,id.vars = 'IID',na.rm = F)
  colnames(df.sub.3) <- c('IID','Tissue','Top')
  df.sub.4 <- dcast(df.sub.3,IID~Tissue)
  rownames(df.sub.4) <- df.sub.4[,1]; df.sub.4 <- df.sub.4[,-1]
  
  tis.type <- colnames(df.sub.4)
  for (i in 1:(length(tis.type)-1)) {
    tis1 <- tis.type[i]
    for (j in (i+1):length(tis.type)) {
      tis2 <- tis.type[j]
      # x <- subset(df.sub.4,!is.na(df.sub.4[,tis1]) & !is.na(`Artery - Coronary`))[,c('Breast - Mammary Tissue','Artery - Coronary')]
      # y <- nrow(subset(x,x[,1]==1 & x[,2]==1))
      # z <- nrow(subset(x,x[,1]==1 | x[,2]==1))
      x <- subset(df.sub.4,!is.na(df.sub.4[,tis1]) & !is.na(df.sub.4[,tis2]))[,c(tis1,tis2)]
      y <- nrow(subset(x,x[,tis1]==1 & x[,tis2]==1))
      z <- nrow(subset(x,x[,tis1]==1 | x[,tis2]==1))
      p <- tryCatch(fisher.test(x[,1],x[,2])$p.value, error=function(e) {NA})
      df.tmp <- data.frame(Tis1=tis1,Tis2=tis2,Prop=y/z,N=z,Pval=p)
      if (i==1 & j==2) {
        df.save <- df.tmp
      } else {
        df.save <- rbind(df.save,df.tmp)
      }
    }
  }
  df.save$Cell <- CellType.unique[k]
  df.save.list[[k]] <- df.save
}
df.save <- do.call(rbind,df.save.list)
df.save$Pval.fdr <- p.adjust(df.save$Pval,method = 'fdr')
subset(df.save,Pval < 0.05/nrow(df.save))
subset(df.save, Pval.fdr < 0.05)
df.save[order(df.save$Prop,decreasing = T),][1:10,]
df.save[order(df.save$Pval,decreasing = F),][1:10,]


subset(df.save,Cell=='CD8Sum' & (Tis1=='Whole Blood' | Tis2=='Whole Blood'))


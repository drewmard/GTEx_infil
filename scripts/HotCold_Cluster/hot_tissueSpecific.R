# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: assess "hot" clusters across tissues within individuals


library(ggplot2)
library(data.table)
library(pheatmap)
library(matrixStats)
library(cowplot)

df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/ConsensusHotColdAssignments.txt',data.table = F,stringsAsFactors = F)
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
  df.sub <- subset(df,CellType==CellType.unique[i])
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
    labs(x='Number of "hot" tissue samples',y='Number of individuals',title=paste0(CellType.names[i],' (N=',length(unique(df.sub$Tissue)),'; n=',val2[[i]],')'))+ 
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5)) +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
}

plot_grid(g[[1]],g[[2]],g[[3]],g[[4]],
          g[[5]],g[[6]],g[[7]],nrow=2)
mean(as.numeric(val))
as.numeric(lapply(dataf,function(x) {median(x$N)}))


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
  x <- apply(df.sub.5,2,sum); df.sub.5 <- df.sub.5[,-which(x==0)]
  f <- paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/HotCold_Cluster/TissueHotHeatMap.',CellType.unique[i],'.png')
  x <-3; p <- pheatmap(df.sub.5,cluster_cols = F,cluster_rows=T,show_colnames = F,
                       legend=F,treeheight_row = 50, treeheight_col = 0,border_color='black',
                       cellwidth = rel(1.2*x),cellheight=rel(4*x),filename=f)
}

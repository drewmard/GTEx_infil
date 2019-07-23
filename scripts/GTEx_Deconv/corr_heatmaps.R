# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script for pairwise heatmaps of GTEx data (Supp fig 2)

library(data.table)
library(ggplot2)

# required files:
df.xcell <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table=F,stringsAsFactors = F)
df.ciber.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.ciber.abs <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

# initialize
dataf <- data.frame(ID=unique(df.ciber.abs$ID))
s <- data.frame(table(df.ciber.rel$SMTSD))

# merge into large dataframe, where each row is new individual and each column is infiltration phenotype
# 212 phenotypes: 53 cell types in 4 tissues
for (i in 1:length(s[,1])) {
  tis <- s[i,1]
  df.sub <- subset(df.ciber.rel,SMTSD %in% tis)[,c('ID','T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
  colnames(df.sub)[-1] <- paste(tis,c('CD8','CD4','Macrophage','Neutrophil'))
  dataf <- merge(dataf,df.sub,by='ID',all=T)
}

# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
g <- ggplot(cor.melt,aes(Var1,Var2,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
                       name="Correlation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title=element_text(hjust=0.5)) + 
  coord_fixed() +
  labs(title='CIBERSORT - Relative')
g

cor.mat.cb_rel <- cor.mat

# TO SAVE PLOTS:
# scale <- 3
# tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','deconv_corr_heatmap_cibersort_rel','.png'),width=700*scale,height=700*scale,res=100*scale,units="px")
# print(g)
# dev.off()

#########################

dataf <- data.frame(ID=unique(df.ciber.abs$ID))
for (i in 1:length(s[,1])) {
  tis <- s[i,1]
  df.sub <- subset(df.ciber.abs,SMTSD %in% tis)[,c('ID','T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
  colnames(df.sub)[-1] <- paste(tis,c('CD8','CD4','Macrophage','Neutrophil'))
  dataf <- merge(dataf,df.sub,by='ID',all=T)
}
# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
g<-ggplot(cor.melt,aes(Var1,Var2,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
                       name="Correlation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title=element_text(hjust=0.5)) + 
  coord_fixed() +
  labs(title='CIBERSORT - Absolute')

cor.mat.cb_abs <- cor.mat

# TO SAVE PLOTS:
# scale <- 3
# tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','deconv_corr_heatmap_cibersort_abs','.png'),width=700*scale,height=700*scale,res=100*scale,units="px")
# print(g)
# dev.off()


########################################

dataf <- data.frame(ID=unique(df.ciber.abs$ID))
for (i in 1:length(s[,1])) {
  tis <- s[i,1]
  df.sub <- subset(df.xcell,SMTSD %in% tis)[,c('ID','CD8Sum','CD4Sum','MacrophageSum','Neutrophils')]
  colnames(df.sub)[-1] <- paste(tis,c('CD8','CD4','Macrophage','Neutrophil'))
  dataf <- merge(dataf,df.sub,by='ID',all=T)
}
# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
g<-ggplot(cor.melt,aes(Var1,Var2,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
                       name="Correlation") +
  theme_minimal()+ # minimal theme
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title=element_text(hjust=0.5)) + 
  coord_fixed() +
  labs(title='xCell')

cor.mat.xcell <- cor.mat

# TO SAVE PLOTS:
# scale <- 3
# tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','deconv_corr_heatmap_xcell','.png'),width=700*scale,height=700*scale,res=100*scale,units="px")
# print(g)
# dev.off()

######################################################################

# mean correlation values:
mean(abs(cor.mat.cb_abs),na.rm=T)
mean(abs(cor.mat.cb_rel),na.rm=T)
mean(abs(cor.mat.xcell),na.rm=T)
mean((cor.mat.cb_abs),na.rm=T)
mean((cor.mat.cb_rel),na.rm=T)
mean((cor.mat.xcell),na.rm=T)

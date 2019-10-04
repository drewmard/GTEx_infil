# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script for pairwise heatmaps of GTEx data (Supp fig 3)

library(data.table)
library(ggplot2)

# required files:
df.xcell <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table=F,stringsAsFactors = F)
df.ciber.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.ciber.abs <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

# initialize
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


# infiltration phenotypes
infil_pheno <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_phenotypes.txt',data.table = F,stringsAsFactors = F)
infil_pheno <- infil_pheno[order(infil_pheno$tissue),]
# infil_pheno$phenotype <- paste(infil_pheno$tissue,infil_pheno$cell,sep = '-')

# merge into large dataframe, where each row is new individual and each column is infiltration phenotype
dataf <- data.frame(ID=unique(df.ciber.abs$ID))
for (i in 1:nrow(infil_pheno)) {
  # tis <- s[i,1]
  tis <- infil_pheno$tissue[i]
  df.sub <- subset(df.ciber.rel,SMTSD %in% tis)[,c('ID',infil_pheno$cell[i])]
  colnames(df.sub)[-1] <- paste(tis,infil_pheno$cell[i])
  dataf <- merge(dataf,df.sub,by='ID',all=T)
}

# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
cor.mat.cb_rel <- cor.mat
cor.mat.melt.cb_rel <- cor.melt
g.cb_rel <- ggplot(cor.mat.melt.cb_rel,aes(Var1,Var2,fill=value)) +
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
  guides(fill=F) +
  coord_fixed() +
  labs(title='CIBERSORT - Relative')
g.cb_rel


#########################

dataf <- data.frame(ID=unique(df.ciber.abs$ID))
for (i in 1:nrow(infil_pheno)) {
  # tis <- s[i,1]
  tis <- infil_pheno$tissue[i]
  df.sub <- subset(df.ciber.abs,SMTSD %in% tis)[,c('ID',infil_pheno$cell[i])]
  colnames(df.sub)[-1] <- paste(tis,infil_pheno$cell[i])
  dataf <- merge(dataf,df.sub,by='ID',all=T)
}
# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
cor.mat.cb_abs <- cor.mat
cor.mat.melt.cb_abs <- cor.melt
g.cb_abs<-ggplot(cor.mat.melt.cb_abs,aes(Var1,Var2,fill=value)) +
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
  guides(fill=F) +
  coord_fixed() +
  labs(title='CIBERSORT - Absolute')
g.cb_abs


########################################

dataf <- data.frame(ID=unique(df.ciber.abs$ID))
for (i in 1:nrow(infil_pheno)) {
  tis <- infil_pheno$tissue[i]
  df.sub <- subset(df.xcell,SMTSD %in% tis)[,c('IID',subset(cellTypes.df,ciber==infil_pheno$cell[i])$xcell)]
  colnames(df.sub)[-1] <- paste(tis,infil_pheno$cell[i])
  dataf <- merge(dataf,df.sub,by.x='ID',by.y='IID',all=T)
}
# pairwise correlation
cor.mat <- cor(dataf[,-1],use='p')

# plot:
cor.melt <- melt(cor.mat)
cor.mat.xcell <- cor.mat
cor.mat.melt.xcell <- cor.melt
g.xcell<-ggplot(cor.mat.melt.xcell,aes(Var1,Var2,fill=value)) +
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
  guides(fill=F) +
  coord_fixed() +
  labs(title='xCell');g.xcell

# library(cowplot)
# plot_grid(g.cb_rel,g.cb_abs,g.xcell,ncol=3)
# 
######################################################################

res <- cor.test(cor.mat.melt.cb_rel$value,cor.mat.melt.cb_abs$value,use='p'); c(res$estimate,res$p.value)

df <- data.frame(cb_rel=cor.mat.melt.cb_rel$value,cb_abs=cor.mat.melt.cb_abs$value,xcell=cor.mat.melt.xcell$value)

g1 <- ggplot(df,aes(x=cb_rel,y=cb_abs)) + geom_point() + geom_smooth(method='lm',se=F,col='red') + labs(x='CIB-Rel',y='CIB-Abs') + theme_bw() + theme(panel.grid=element_blank())
g2 <- ggplot(df,aes(x=xcell,y=cb_abs)) + geom_point() + geom_smooth(method='lm',se=F,col='red') + labs(x='xCell',y='CIB-Abs') + theme_bw() + theme(panel.grid=element_blank())
g3 <- ggplot(df,aes(x=xcell,y=cb_rel)) + geom_point() + geom_smooth(method='lm',se=F,col='red') + labs(x='xCell',y='CIB-Rel') + theme_bw()+ theme(panel.grid=element_blank())
plot_grid(g1,g2,g3,ncol=3)

plot(cor.mat.melt.cb_rel$value,cor.mat.melt.cb_abs$value)

res <- cor.test(cor.mat.melt.xcell$value,cor.mat.melt.cb_abs$value,use='p'); c(res$estimate,res$p.value)
res <- cor.test(cor.mat.melt.xcell$value,cor.mat.melt.cb_rel$value,use='p'); c(res$estimate,res$p.value)

# i <- which(cor.mat.melt.cb_rel$value > 0.5)
# i <- which(cor.mat.melt.cb_abs$value > 0.2)
# much stronger correlation between xcell & cibersort when considering highly correlation infiltration patterns
i <- which(cor.mat.melt.xcell$value > 0.2)
res <- cor.test(cor.mat.melt.cb_rel$value[i],cor.mat.melt.cb_abs$value[i],use='p'); c(res$estimate,res$p.value)
res <- cor.test(cor.mat.melt.xcell$value[i],cor.mat.melt.cb_abs$value[i],use='p'); c(res$estimate,res$p.value)
res <- cor.test(cor.mat.melt.xcell$value[i],cor.mat.melt.cb_rel$value[i],use='p'); c(res$estimate,res$p.value)

# mean correlation values:
mean(abs(cor.mat.cb_abs),na.rm=T)
mean(abs(cor.mat.cb_rel),na.rm=T)
mean(abs(cor.mat.xcell),na.rm=T)
mean((cor.mat.cb_abs),na.rm=T)
mean((cor.mat.cb_rel),na.rm=T)
mean((cor.mat.xcell),na.rm=T)

library(cowplot)
plot_grid(g.cb_rel,g.cb_abs,g.xcell,ncol=3)



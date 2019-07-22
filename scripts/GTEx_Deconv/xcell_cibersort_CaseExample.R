# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script for specific case examples presented in supp fig 3

library(data.table)
library(stringr)
library(ggplot2)
library(cowplot)

# files
df.xcell <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table=F,stringsAsFactors = F)
df.infil.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.infil.abs <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

df = data.frame(abs=subset(df.infil.abs,SMTSD=='Lung')$Neutrophils,
                rel=subset(df.infil.rel,SMTSD=='Lung')$Neutrophils,
                xcell=subset(df.xcell,SMTSD=='Lung')$Neutrophils
)

# Individuals 1,2

x1 = subset(df.infil.rel,SMTSD=='Lung')[c(159,195),c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
x2 = subset(df.infil.abs,SMTSD=='Lung')[c(159,195),c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
x3 = subset(df.xcell,SMTSD=='Lung')[c(159,195),c('CD8Sum','CD4Sum','MacrophageSum','Neutrophils')]

res = data.frame(Method=rep(c('CIBERSORT-Rel','CIBERSORT-Abs','xCell'),each=2),Indiv=rep(c('GTEX-14LLW','GTEX-17F96'),3),Neutrophils=c(x1$Neutrophils,x2$Neutrophils,x3$Neutrophils))
res$Method <- as.factor(res$Method)
res$Indiv <- as.factor(res$Indiv)

g1 <- ggplot(res[1:2,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='CIBERSORT-Relative') +
  scale_fill_manual(values=c('steelblue2','orange2'))
g2 <- ggplot(res[3:4,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='CIBERSORT-Absolute') +
  scale_fill_manual(values=c('steelblue2','orange2'))
g3 <- ggplot(res[5:6,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='xCell') +
  scale_fill_manual(values=c('steelblue2','orange2'))

plot_grid(g1,g2,g3,ncol=3)

# Individuals 3,4

x1 = subset(df.infil.rel,SMTSD=='Lung')[c(269,405),c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
x2 = subset(df.infil.abs,SMTSD=='Lung')[c(269,405),c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')]
x3 = subset(df.xcell,SMTSD=='Lung')[c(269,405),c('CD8Sum','CD4Sum','MacrophageSum','Neutrophils')]

res = data.frame(Method=rep(c('CIBERSORT-Rel','CIBERSORT-Abs','xCell'),each=2),Indiv=rep(c('GTEX-14LLW','GTEX-17F96'),3),Neutrophils=c(x1$Neutrophils,x2$Neutrophils,x3$Neutrophils))
res$Method <- as.factor(res$Method)
res$Indiv <- as.factor(res$Indiv)

g1 <- ggplot(res[1:2,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='CIBERSORT-Relative') +
  scale_fill_manual(values=c('steelblue2','orange2'))
g2 <- ggplot(res[3:4,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='CIBERSORT-Absolute') +
  scale_fill_manual(values=c('steelblue2','orange2'))
g3 <- ggplot(res[5:6,],aes(x=Method,fill=Indiv,y=Neutrophils)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(legend.position = 'none',axis.title.x = element_blank()) +
  scale_x_discrete(labels=NULL) +
  labs(title='xCell') +
  scale_fill_manual(values=c('steelblue2','orange2'))
plot_grid(g1,g2,g3,ncol=3)

# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: perform t-sne clustering on different tissues to see if any clusters develop

library(data.table)
library(stringr)
library(tsne)
library(ggplot2)

deconv_data.Cibersort <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
# aggregate(deconv_data.Cibersort$`B cells naive`,list(deconv_data.Cibersort$SMTSD),length)
Tissues <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/TISSUES_n70.txt',data.table = F,stringsAsFactors = F,sep=',')[,1]
x <- subset(deconv_data.Cibersort,SMTSD==Tissues[5])
x1 <- (x)[,3:24]
tsne <- tsne(x1)
tsne1 <- data.frame('Component.1'=tsne[,1],'Component.2'=tsne[,2])
i = 4; colnames(x)[i]
g1 <- ggplot(tsne1,aes(x=`Component.1`,y=`Component.2`,col=x1[,i])) + geom_point() + ggtitle(Tissues[5]) +
  labs(col=colnames(x)[i])+ theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  scale_color_continuous(low='lightblue',high='red')

deconv_data.Cibersort <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
Tissues <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/TISSUES_n70.txt',data.table = F,stringsAsFactors = F,sep=',')[,1]
x <- subset(deconv_data.Cibersort,SMTSD==Tissues[23])
x2 <- (x)[,3:24]
tsne <- tsne(x2)
tsne2 <- data.frame('Component.1'=tsne[,1],'Component.2'=tsne[,2])
i = 4; colnames(x)[i]
g2 <- ggplot(tsne2,aes(x=`Component.1`,y=`Component.2`,col=x2[,i])) + geom_point() + ggtitle(Tissues[23]) +
  labs(col=colnames(x)[i]) + theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  scale_color_continuous(low='lightblue',high='red')

library(data.table)
library(stringr)
deconv_data.Cibersort <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
Tissues <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/TISSUES_n70.txt',data.table = F,stringsAsFactors = F,sep=',')[,1]
x <- subset(deconv_data.Cibersort,SMTSD==Tissues[9])
x3 <- (x)[,3:24]
tsne <- tsne(x3)
tsne3 <- data.frame('Component.1'=tsne[,1],'Component.2'=tsne[,2])
i = 4; colnames(x)[i]
g3 <- ggplot(tsne3,aes(x=`Component.1`,y=`Component.2`,col=x3[,i])) + geom_point() + ggtitle(Tissues[9]) +
  labs(col=colnames(x)[i])+ theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  scale_color_continuous(low='lightblue',high='red')

library(data.table)
library(stringr)
deconv_data.Cibersort <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
Tissues <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/TISSUES_n70.txt',data.table = F,stringsAsFactors = F,sep=',')[,1]
x <- subset(deconv_data.Cibersort,SMTSD==Tissues[32])
x4 <- (x)[,3:24]
tsne <- tsne(x4)
tsne4 <- data.frame('Component.1'=tsne[,1],'Component.2'=tsne[,2])
i = 4; colnames(x)[i]
g4 <- ggplot(tsne4,aes(x=`Component.1`,y=`Component.2`,col=x4[,i])) + geom_point() + ggtitle(Tissues[32]) +
  labs(col=colnames(x)[i]) + theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
  scale_color_continuous(low='lightblue',high='red')

library(cowplot)
plot_grid(g1,g2,g3,g4,ncol=2)


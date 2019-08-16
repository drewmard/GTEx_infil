# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: create venn diagrams displaying joint (EBM) analysis versus separate analyses

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
library(limma)
library(data.table)

# all analyses: p < 5e-8
# f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.sig2.txt'
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/GWAS_iQTL.txt'
df.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
df.sub.save$Pval_Brown_sig <- as.numeric(df.sub.save$Pval_Brown < 5e-8)
df.sub.save$Pval_Brown_Abs_sig <- as.numeric(df.sub.save$Pval_Brown_Abs < 5e-8)
df.sub.save$p_lrt_cibersort_rel_sig <- as.numeric(df.sub.save$p_lrt_cibersort_rel < 5e-8)
df.sub.save$p_lrt_cibersort_abs_sig <- as.numeric(df.sub.save$p_lrt_cibersort_abs < 5e-8)
df.sub.save$p_lrt_cibersort_xCell_sig <- as.numeric(df.sub.save$p_lrt_xCell < 5e-8)
x <- aggregate(df.sub.save[,c(15,17:19)],by=list(df.sub.save$tissue,df.sub.save$cell),max)
x <- subset(x,!(Pval_Brown_sig==0 & p_lrt_cibersort_rel_sig==0 &
              p_lrt_cibersort_abs_sig==0 & p_lrt_cibersort_xCell_sig==0))
a1 <- vennCounts(x[,-(1:2)])
colnames(a1)[1:4] <- c('EBM','Cibersort-Rel','Cibersort-Abs','xCell')
vennDiagram(a1,circle.col=c('steelblue1','steelblue2','steelblue3','steelblue4'))

# EBM p < 5e-8, separate: p < 1e-5
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/GWAS_iQTL.txt'
df.sub.save <- fread(f,data.table = F,stringsAsFactors = F)
df.sub.save$Pval_Brown_sig <- as.numeric(df.sub.save$Pval_Brown < 5e-8)
df.sub.save$Pval_Brown_Abs_sig <- as.numeric(df.sub.save$Pval_Brown_Abs < 5e-8)
df.sub.save$p_lrt_cibersort_rel_sig <- as.numeric(df.sub.save$p_lrt_cibersort_rel < 1e-5)
df.sub.save$p_lrt_cibersort_abs_sig <- as.numeric(df.sub.save$p_lrt_cibersort_abs < 1e-5)
df.sub.save$p_lrt_cibersort_xCell_sig <- as.numeric(df.sub.save$p_lrt_xCell < 1e-5)
x <- aggregate(df.sub.save[,c(15,17:19)],by=list(df.sub.save$tissue,df.sub.save$cell),max)
x <- subset(x,!(Pval_Brown_sig==0 & p_lrt_cibersort_rel_sig==0 &
                  p_lrt_cibersort_abs_sig==0 & p_lrt_cibersort_xCell_sig==0))
a2 <- vennCounts(x[,-(1:2)])
colnames(a2)[1:4] <- c('EBM','Cibersort-Rel','Cibersort-Abs','xCell')
vennDiagram(a2,circle.col=c('steelblue1','steelblue2','steelblue3','steelblue4'))

# age & sex venn diagram
df.coef <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/AgeSex_results.txt',data.table = F,stringsAsFactors = F)
df.coef$xcell.fdr <- as.numeric(p.adjust(df.coef$p_age_xcell,method='fdr')<0.1 | p.adjust(df.coef$p_sex_xcell,method='fdr')<0.1)
df.coef$cib.abs.fdr <- as.numeric(p.adjust(df.coef$p_age_cib.abs,method='fdr')<0.1 | p.adjust(df.coef$p_sex_cib.abs,method='fdr')<0.1)
df.coef$cib.rel.fdr <- as.numeric(p.adjust(df.coef$p_age_cib.rel,method='fdr')<0.1 | p.adjust(df.coef$p_sex_cib.rel,method='fdr')<0.1)
df.coef$brown.fdr <- as.numeric(df.coef$p_age_brown.fdr<0.1 | df.coef$p_sex_brown.fdr<0.1)
df.coef$xcell.fdr[is.na(df.coef$xcell.fdr)] <- 0
df.coef$cib.abs.fdr[is.na(df.coef$cib.abs.fdr)] <- 0
df.coef$cib.rel.fdr[is.na(df.coef$cib.rel.fdr)] <- 0
df.coef$brown.fdr[is.na(df.coef$brown.fdr)] <- 0
x <- df.coef
a <- vennCounts(x[,c('brown.fdr','cib.rel.fdr','cib.abs.fdr','xcell.fdr')])
colnames(a)[1:4] <- c('EBM','Cibersort-Rel','Cibersort-Abs','xCell')
vennDiagram(a,circle.col=c('steelblue1','steelblue2','steelblue3','steelblue4'))

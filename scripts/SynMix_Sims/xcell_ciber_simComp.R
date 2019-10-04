# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for comparing xcell vs cibersort in simulation performance

library(data.table)
ciber.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/SynMix_1.CIBERSORT.ABS-F.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
ciber.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/SynMix_1.CIBERSORT.ABS-T.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput_merged.txt',data.table = F,stringsAsFactors = F)

xcell.normalized=FALSE
if (xcell.normalized) {xcell[,2:65] <- xcell[,2:65]/apply(xcell[,2:65],1,sum)}
xcell$CD4_Tcells_Sum <- apply(xcell[,c('CD4+ memory T-cells',
                                       'CD4+ naive T-cells',
                                       'CD4+ T-cells',
                                       'CD4+ Tcm',
                                       'CD4+ Tem',
                                       'Th1 cells',
                                       'Th2 cells',
                                       'Tregs'
                                       )],1,sum)
xcell$CD8_Tcells_Sum <- apply(xcell[,c('CD8+ naive T-cells',
                                       'CD8+ T-cells',
                                       'CD8+ Tcm',
                                       'CD8+ Tem')],1,sum)
ciber.rel$CD4_Tcells_Sum <- apply(ciber.rel[,c('T cells CD4 naive',
                                           'T cells CD4 memory resting',
                                           'T cells CD4 memory activated',
                                           'T cells follicular helper',
                                           'T cells regulatory (Tregs)')],1,sum)
ciber.abs$CD4_Tcells_Sum <- apply(ciber.abs[,c('T cells CD4 naive',
                                           'T cells CD4 memory resting',
                                           'T cells CD4 memory activated',
                                           'T cells follicular helper',
                                           'T cells regulatory (Tregs)')],1,sum)
ciber.abs$CD8_Tcells_Sum <- ciber.abs[,'T cells CD8']
ciber.rel$CD8_Tcells_Sum <- ciber.rel[,'T cells CD8']

ciber.abs$`CD8_Tcells_Sum` <- ciber.abs$`T cells CD8`
ciber.rel$`CD8_Tcells_Sum` <- ciber.rel$`T cells CD8`

# xcell[,2:65] <- xcell[,2:65]/apply(xcell[,2:65],1,sum)
# xcell$CD4_Tcells_Sum <- apply(xcell[,7:11],1,sum)
# xcell$CD8_Tcells_Sum <- apply(xcell[,12:15],1,sum)
# 
# ciber.abs$`CD8_Tcells_Sum` <- ciber.abs$`T cells CD8`
# ciber.rel$`CD8_Tcells_Sum` <- ciber.rel$`T cells CD8`

cor(xcell$`CD4_Tcells_Sum`,xcell$CD4_Naive_Actual_Pct)
cor(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual_Pct)
cor(xcell$`CD4_Tcells_Sum`,xcell$CD4_Naive_Actual)
cor(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual)

cor(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual_Pct)
cor(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual_Pct)
cor(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual)
cor(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual)

cor(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual_Pct)
cor(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual_Pct)
cor(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual)
cor(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual)

cor(ciber.abs$`CD8_Tcells_Sum`,xcell$`CD8_Tcells_Sum`)
cor(ciber.abs$`CD4_Tcells_Sum`,xcell$`CD4_Tcells_Sum`)

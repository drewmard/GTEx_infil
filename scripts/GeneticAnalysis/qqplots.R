# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to draw qq plots

library(data.table); library(ggplot2); library(qqman)

# Arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) # what pheno to look at?

# load data
workdir <- '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/'
f <- paste0(workdir,'GTEx.pheno',i,'.ALL_EBM.txt')

pvector <- df$Pval_Brown
n <- length(pvector)+1
exp.x <- -log10((rank(pvector, ties.method="first")-.5)/n)
pvalues <- -log10(pvector)
g <- ggplot(data.frame(exp=exp.x,obs=pvalues),aes(x=exp,y=obs)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
  theme_bw() + theme(panel.grid=element_blank(),plot.title = element_text(hjust=0.5)) + 
  scale_y_continuous(breaks=seq(0,max(pvalues),by=2)) +
  scale_x_continuous(breaks=seq(0,max(exp.x),by=2)) +
  labs(x = expression(Expected ~ ~-log[10](italic(p))), y = expression(Observed ~ ~-log[10](italic(p)))) +
  ggtitle('Immune cell type in tissue');# g

# saving:
# x <- 3
# tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','Skin_-_Sun_Exposed_(Lower_leg).1','.png'),width=400*x,height=400*x,res=100*x,units="px")
# print(g)
# dev.off()

library(data.table)
library(ggplot2)
library(cowplot)

i=1
id[[1]] <- 166
workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/'
f <- paste0(workdir,'output/GeneticAnalysis/GWAS/','GTEx.pheno',id[[i]],'.ALL_EBM.txt')
df <- fread(f,data.table = F,stringsAsFactors = F)

df.qq <- list(); g.qq <- list()
col_names <- c('p_lrt_cibersort_rel',
               "p_lrt_cibersort_abs",
               "p_lrt_xCell")
for (j in 1:3) {
  # pvector <- df[1:10000,col_names[j]]
  pvector <- df[,col_names[j]]
  n <- length(pvector)+1
  exp.x <- -log10((rank(pvector, ties.method="first")-.5)/n)
  pvalues <- -log10(pvector)
  df.qq[[j]] <- data.frame(exp=exp.x,obs=pvalues)
  g.qq[[j]] <- ggplot(df.qq[[j]],aes(x=exp,y=obs)) + geom_point() + geom_abline(slope=1,intercept=0,col='red') + 
    theme_bw() + theme(panel.grid=element_blank(),plot.title = element_text(hjust=0.5)) + 
    scale_y_continuous(breaks=seq(0,max(pvalues),by=2)) +
    scale_x_continuous(breaks=seq(0,max(exp.x),by=2)) +
    labs(x = expression(Expected ~ ~-log[10](italic(p))), y = expression(Observed ~ ~-log[10](italic(p)))) #+
  
}

figure_panel <- plot_grid(g.qq[[1]],g.qq[[2]],g.qq[[3]],ncol=3)

# tiff('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/figure_panel.png',width=4500,height=5000,res=400,units="px")
tiff('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/qq.png',width=4500,height=1600,res=400,units="px")
print(figure_panel)
dev.off()

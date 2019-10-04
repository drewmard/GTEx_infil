# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for plotting ieQTL/eQTL enrichment results

library(data.table)
library(ggplot2)

# 1 sided test

df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/eQTL_enrich2_1.txt',data.table=F,stringsAsFactors = F)
df$Analysis <- paste0(gsub('_',' ',df$tissue),': ',df$cell)
df[order(df$p,decreasing = T),]
df$eqtl.enriched <- as.numeric(df$eQTL_obs_ct/df$N > df$mean.p)
df$sig <- 1
df$sig[df$p > 0.05/nrow(df)] <- 0
g.matched <- ggplot(df,aes(x=reorder(Analysis,-p),y=-log10(p),col=as.factor(sig))) + 
  geom_point(size=rel(0.7)) +
  geom_hline(yintercept=-log10(0.05/nrow(df)),col='red',lty='dashed') +
  scale_color_manual(values = c('black','brown2'),labels=c('No significance','Significant enrichment')) +
  labs(x='Infiltration Phenotype',y=expression(~-log[10](italic(p)))) +
  ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
  theme_bw() + theme(panel.grid = element_blank(),
                     # axis.text.x = element_text(angle=60,hjust=1,size=rel(0.6)),
                     axis.text.x = element_blank(),
                     legend.title = element_blank(),legend.background = element_blank(),
                     legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
  theme(legend.position = c(0.1,0.85)); g.matched
x <- 3
tiff('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/eQTL_enrichment.matched.1.png',width = 1000*x,height=250*x,res=100*x)
print(g.matched)
dev.off()


# 2 sided test

df$eqtl.enriched2 <- df$eqtl.enriched + 1
df$eqtl.enriched2[df$p2 > 0.05/nrow(df)] <- 0
g.matched <- ggplot(df,aes(x=reorder(Analysis,-p2),y=-log10(p2),col=as.factor(eqtl.enriched2))) + 
  geom_point(size = rel(0.7)) +
  geom_hline(yintercept=-log10(0.05/nrow(df)),col='red',lty='dashed') +
  scale_color_manual(values = c('black','skyblue1','brown2'),labels=c('Non-significant','Fewer eQTLs than expected','Greater eQTLs than expected')) +
  labs(x='Infiltration phenotype',y=expression(~-log[10](italic(p)))) +
  ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
  theme_bw() + theme(panel.grid = element_blank(),
                     # axis.text.x = element_text(angle=60,hjust=1,size=rel(0.6)),
                     axis.text.x = element_blank(),
                     legend.title = element_blank(),legend.background = element_blank(),
                     legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
  theme(legend.position = c(0.1,0.85)); g.matched

x <- 3
tiff('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/eQTL_enrichment.matched.2.png',width = 1000*x,height=250*x,res=100*x)
print(g.matched)
dev.off()



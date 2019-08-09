library(data.table)
library(ggplot2)

# df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/020819_eQTL_enrichment.1.csv',data.table=F,stringsAsFactors = F)
df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis2/eQTL_enrich2_1.txt',data.table=F,stringsAsFactors = F)
# df$cell[df$cell=='MacrophageSum'] <- 'Macrophages'
# df$cell[df$cell=='T cells CD8'] <- 'CD8+ T cells'
# df$cell[df$cell=='CD4_Tcells'] <- 'CD4+ T cells'
# df$tissue2 <- gsub('_',' ',df$tissue)
df$Analysis <- paste0(gsub('_',' ',df$tissue),': ',df$cell)
g.matched <- ggplot(df,aes(x=reorder(Analysis,-p),y=-log10(p))) + 
  geom_point() +
  geom_hline(yintercept=-log10(0.05/nrow(df)),col='red',lty='dashed') +
  scale_color_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +
  labs(x='Tissue and Immune Cell Type',y=expression(~-log[10](italic(p)))) +
  ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
  # ggtitle(expression(Observed ~ ~-log[10](italic(p)))) +
  # scale_color_discrete(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle=80,hjust=1,size=rel(0.8)),
                     legend.title = element_blank(),legend.background = element_blank(),
                     legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
  theme(legend.position = c(0.1,0.85)); g.matched

mean(df$p < 0.05/nrow(df),na.rm=T)


df2 <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/020819_eQTL_enrichment.1.csv',data.table=F,stringsAsFactors = F)
df2$Cell[df2$Cell=='MacrophageSum'] <- 'Macrophages'
df2$Cell[df2$Cell=='T cells CD8'] <- 'CD8+ T cells'
df2$Cell[df2$Cell=='CD4_Tcells'] <- 'CD4+ T cells'
df2$tissue2 <- gsub('_',' ',df2$Tissue)
df2$Analysis <- paste0(gsub('_',' ',df2$tissue2),': (',df2$Cell,')')
df <- merge(df,df2[,c('Analysis','p.value','eQTL.enriched')],by='Analysis')

g.chisq <- ggplot(df,aes(x=reorder(Analysis,-p.value),y=-log10(p.value),col=as.factor(eQTL.enriched))) + 
  geom_point() +
  scale_color_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +
  labs(x='Tissue and Immune Cell Type',y=expression(~-log[10](italic(p)))) +
  ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
  # ggtitle(expression(Observed ~ ~-log[10](italic(p)))) +
  # scale_color_discrete(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle=80,hjust=1,size=rel(0.8)),
                     legend.title = element_blank(),legend.background = element_blank(),
                     legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
  theme(legend.position = c(0.1,0.85)); g

g.matched <- ggplot(df,aes(x=reorder(Analysis,-pval),y=-log10(pval),col=as.factor(eQTL.enriched))) + 
  geom_point() +
  scale_color_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +
  labs(x='Tissue and Immune Cell Type',y=expression(~-log[10](italic(p)))) +
  ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
  # ggtitle(expression(Observed ~ ~-log[10](italic(p)))) +
  # scale_color_discrete(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle=80,hjust=1,size=rel(0.8)),
                     legend.title = element_blank(),legend.background = element_blank(),
                     legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
  theme(legend.position = c(0.1,0.85)); g

# g <- ggplot(df,aes(x=reorder(Analysis,-p.value),y=-log10(p.value),col=as.factor(eQTL.enriched))) + 
#   geom_point() +
#   geom_point(aes(y=-log10(pval)),shape=2) +
#   scale_color_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +
#   labs(x='Tissue and Immune Cell Type',y=expression(~-log[10](italic(p)))) +
#   ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
#   # ggtitle(expression(Observed ~ ~-log[10](italic(p)))) +
#   # scale_color_discrete(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),) +
#   theme_bw() + theme(panel.grid = element_blank(),
#                      axis.text.x = element_text(angle=80,hjust=1,size=rel(0.8)),
#                      legend.title = element_blank(),legend.background = element_blank(),
#                      legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
#   theme(legend.position = c(0.1,0.85)); g
# g <- ggplot(df,aes(x=reorder(Analysis,-pval),y=-log10(pval),col=as.factor(eQTL.enriched))) + 
#   geom_point() +
#   geom_point(aes(y=-log10(p.value)),shape=2) +
#   # scale_color_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +
#   scale_colour_manual(values = c('skyblue1','brown2'),labels=c('Fewer eQTLs than expected','Greater eQTLs than expected')) +   
#   scale_shape_manual(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),
#                      values = c(1, 2)) +
#   labs(x='Tissue and Immune Cell Type',y=expression(~-log[10](italic(p)))) +
#   ggtitle(bquote('Enrichment of tissue-specific eQTLs in top GWAS hits ( p < '~10^-5~')')) +
#   # ggtitle(expression(Observed ~ ~-log[10](italic(p)))) +
#   # scale_color_discrete(labels=c('Fewer eQTLs than expected','Greater eQTLs than expected'),) +
#   theme_bw() + theme(panel.grid = element_blank(),
#                      axis.text.x = element_text(angle=80,hjust=1,size=rel(0.8)),
#                      legend.title = element_blank(),legend.background = element_blank(),
#                      legend.key = element_blank(),plot.title=element_text(hjust=0.5)) +
#   theme(legend.position = c(0.1,0.85)); g
x <- 3
tiff('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/eQTL_enrichment.matched.png',width = 1000*x,height=480*x,res=100*x)
print(g.matched)
dev.off()
x <- 3
tiff('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/eQTL_enrichment.chisq.png',width = 1000*x,height=480*x,res=100*x)
print(g.chisq)
dev.off()



head(df)
wilcox.test(-log10(subset(df,eQTL.enriched==0)$p.value),-log10(subset(df,eQTL.enriched==1)$p.value),alternative = 'less')
t.test(-log10(subset(df,eQTL.enriched==0)$p.value),-log10(subset(df,eQTL.enriched==1)$p.value),alternative = 'less')

wilcox.test(-log10(subset(df,eQTL.enriched==0)$pval),-log10(subset(df,eQTL.enriched==1)$pval),alternative = 'less')
t.test(-log10(subset(df,eQTL.enriched==0)$pval),-log10(subset(df,eQTL.enriched==1)$pval),alternative = 'less')



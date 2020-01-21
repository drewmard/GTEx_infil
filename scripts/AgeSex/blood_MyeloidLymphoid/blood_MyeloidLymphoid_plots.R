library(data.table)
df.save <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/ratio.df.txt',data.table = F,stringsAsFactors = F)
colnames(df.save)[2:4] <- c('xCell','CIBER','Age')
library(ggplot2)
g1 <- ggplot(df.save,aes(x=CIBER,y=xCell)) + geom_point() + theme_bw() +
  labs(x='Myeloid:Lymphoid Ratio (CIBERSORT; trans resid)',y='Myeloid:Lymphoid Ratio (xCell; trans resid)'); g1
g2 <- ggplot(df.save,aes(x=Age,y=xCell)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1,alpha=0.2,col='orange') + theme_bw() +
  labs(y='Myeloid:Lymphoid Ratio (xCell; trans resid)');
g3 <- ggplot(df.save,aes(x=Age,y=CIBER)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1,alpha=0.2,col='orange') + theme_bw() +
  labs(y='Myeloid:Lymphoid Ratio (CIBERSORT; trans resid)');

library(cowplot)
plot_grid(g1,g2,g3,ncol=3)

cor.test(df.save$xCell,df.save$CIBER,use='p')$p.value
cor.test(as.numeric(as.factor(df.save$Age)),df.save$CIBER,use='p')
cor.test(as.numeric(as.factor(df.save$Age)),df.save$xCell,use='p')

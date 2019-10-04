df <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df <- fread('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.sub <- subset(df,SMTSD=='Breast - Mammary Tissue')
library("tsne")
res <- tsne(df.sub[,3:24])
dataf <- data.frame(tsne1=res[,1],tsne2=res[,2],sex=as.factor(df.sub$SEX))

plot(dataf$tsne1,dataf$tsne2,col=dataf$sex,pch=16,
     xlab='Component 1',ylab='Component 2',main='t-SNE of immune content in male/female breast tissue')

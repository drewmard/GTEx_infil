# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: analysis of genetic influences w/ thyroiditis that are associated w/ commd3 and dnajc1

library(phenoscanner)

# self reported thyroiditis GWAS results from neale lab, downloaded from neale lab ukbiobank website
df <- read.table('/athena/elementolab/scratch/anm2868/GTEx/20002_1428.gwas.imputed_v3.both_sexes.tsv.bgz',stringsAsFactors = F,header=T)

# organize variant name into chr/pos variables
x <- strsplit(df$variant,':')
df$CHR <- unlist(x)[c(T,F,F,F)]
df$POS <- unlist(x)[c(F,T,F,F)]
subset(df,CHR==10 & POS==22337752)

# commd3/dnajc1 promoter/enhancer
region <- c(
  'chr10:22517118-22520043',
  'chr10:22539329-22544453',
  'chr10:22623174-22624900',
  'chr10:22621725-22623130',
  'chr10:22625421-22626070',
  'chr10:22722368-22727938'
)
region.gwas <- data.frame(chromosome=as.numeric(substring(unlist(strsplit(region,':|-'))[c(T,F,F)],4)),
                          start=as.numeric(unlist(strsplit(region,':|-'))[c(F,T,F)]),
                          end=as.numeric(unlist(strsplit(region,':|-'))[c(F,F,T)]),
                          stringsAsFactors = F
)

# query gwas results:
for (i in 1:nrow(region.gwas)) {
  df.sub <- subset(df,as.numeric(CHR)==region.gwas$chromosome[1] & as.numeric(POS) > as.numeric(region.gwas$start[i]) & as.numeric(POS) < as.numeric(region.gwas$end[i]))
  if (i == 1) {
    df.sub.save <- df.sub
  } else {
    df.sub.save <- rbind(df.sub.save,df.sub)
  }
}

# query MAF > 0.01
df.sub.save2 <- subset(df.sub.save,as.numeric(as.character(minor_AF)) > 0.01)

# snps passing bonferroni significance
subset(df.sub.save2,pval < 0.05/nrow(df.sub.save2))

# can also check phenoscanner:
res <- phenoscanner(regionquery = region,pval=0.05,r2=0.8,catalogue = 'GWAS',build=38)
res.sub <- res$results
res.sub[order(as.numeric(res.sub$p),decreasing = F),c('trait','p','n')]

# done
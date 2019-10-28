# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: global enrichment analysis of iQTLs for previous GWAS hits (P < 5e-8)

library(phenoscanner)
library(data.table)

# read in & identify 31 iQTLs w/ p < 5e-8 from 31 infiltration phenotypes w/ genome-wide sig associations
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/Supplementary/SuppTab11.txt'
df <- fread(f,data.table=F,stringsAsFactors = F)
df$pheno <- paste(df$tissue,df$cell,sep=' ')
pheno.uniq <- unique(df$pheno)
SNP.vec <- c()
for (i in 1:length(pheno.uniq)) {
  df.tmp <- subset(df,pheno==pheno.uniq[i])
  i <- which.min(df.tmp$Pval_Brown)[1]
  SNP.vec <- c(SNP.vec,df.tmp$SNP[i])
}

# to map GTEx id to rs id, using downloaded GTEx v6 conversion table
f <- '~/Downloads/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1-22-X_genot_imput_info04_maf01_HWEp1E6_variant_id_lookup.txt'
snp_info <- fread(f,data.table = F,stringsAsFactors = F)
snp_info2 <- subset(snp_info,variant_id %in% SNP.vec)
snp_info2 <- snp_info2[,c('variant_id','rs_id_dbSNP142_GRCh37p13')]
colnames(snp_info2)[2] <- 'rs'
df2 <- merge(df,snp_info2,by.x='SNP',by.y='variant_id')

# previous GWAS hit?
gwas_hits.vec <- c()
for (i in 1:length(pheno.uniq)) {
  n <- nrow(phenoscanner(snpquery=df2$rs[i],r2=0.8,pvalue = 5e-8)$results)
  gwas_hits.vec <- c(gwas_hits.vec,n)
}
sum(gwas_hits.vec > 0)/length(gwas_hits.vec)

# permutations
p.vec <- c()
for (j in 1:10) {
  random_snp <- sample(snp_info$rs_id_dbSNP142_GRCh37p13,length(pheno.uniq),replace = F)
  gwas_hits.vec5 <- c()
  for (i in 1:length(pheno.uniq)) {
    n <- nrow(phenoscanner(snpquery=random_snp[i],r2=0.8,pvalue = 5e-8)$results)
    gwas_hits.vec5 <- c(gwas_hits.vec5,n)
  }
  p.vec <- c(p.vec,sum(gwas_hits.vec5 > 0)/length(gwas_hits.vec5))
  print(p.vec)
}

# p-value enrichment of GWAS hit observed
q <- mean(p.vec)
N <- length(gwas_hits.vec)
x <- sum(gwas_hits.vec > 0)
binom.test(x,N,q,alternative = 'greater')

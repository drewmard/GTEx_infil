library(data.table)
library(stringr)

# LD european population patterns - TAKES A WHILE TO READ IN.
LD <- fread('/athena/elementolab/scratch/anm2868/GTEx/LD_EUR.tsv',data.table=F,stringsAsFactor=F,sep='\t',header=F)
# frq <- fread('/athena/elementolab/scratch/anm2868/GTEx/plink.frq',data.table = F,stringsAsFactors = F)

# need to switch directory
gtex_to_rsid <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_rsid_matched.txt',data.table=F,stringsAsFactors=F,header=F)

# number of snps in LD
LD$LD_SNP_ct <- str_count(LD[,2],';')+1
# LD.old <- LD

# rs to gtex names
# head(LD,1)
# head(gtex_to_rsid,1)
LD <- merge(LD,gtex_to_rsid,by.x='V1',by.y='V2')

# add in freq 
colnames(LD) <- c('rsid','LD_snps','LD_SNP_ct','gtex_snpid')
LD <- merge(LD,frq,by.x='gtex_snpid',by.y='SNP') # switch by.x

# eQTL <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_Analysis_v7_eQTL/',TISSUE.eqtl,'.v7.signif_variant_gene_pairs.txt'),data.table = F)

# match rs id to gtex
# sig_eQTL <- fread('/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/Subset_Pval.eQTL.1.csv',data.table=F,stringsAsFactors = F)
sig_eQTL <- fread('/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/020619_Subset_Pval.eQTL.1.csv',data.table=F,stringsAsFactors = F)
df.to.run = unique(sig_eQTL[,c('tissue','cell','pheno_num')]); rownames(df.to.run) <- 1:nrow(df.to.run)
nperm = 10
pval.vec <- rep(NA,nrow(df.to.run))
for (i in 1:nrow(df.to.run)) {
  print(i)
  ind_lst <- list()
  tis <- df.to.run$tissue[i]; cel <- df.to.run$cell[i]; phenoNum <- df.to.run$pheno_num[i];
  df.tmp = subset(sig_eQTL,tissue==tis & cell==cel & pheno_num==phenoNum)
  df.tmp = merge(df.tmp,LD,by.x='rs',by.y='gtex_snpid')
  df.tmp <- df.tmp[,-which(colnames(df.tmp)=='LD_snps')]
  for (j in 1:nrow(df.tmp)) {
    ind <- which(LD$LD_SNP_ct == (df.tmp$LD_SNP_ct[j]) &
                     LD$MAF > (df.tmp$MAF[j]-0.01) &
                     LD$MAF < (df.tmp$MAF[j]+0.01))
    if (length(ind) <= 1) {
      ind <- which(LD$LD_SNP_ct > (df.tmp$LD_SNP_ct[j]-5) &
                    LD$LD_SNP_ct < (df.tmp$LD_SNP_ct[j]+5) &
                    LD$MAF > (df.tmp$MAF[j]-0.01) &
                    LD$MAF < (df.tmp$MAF[j]+0.01))
    }
    ind <- ind[ind!=match(df.tmp$rs[j],LD$gtex_snpid)]
    ind_lst[[j]] <- sample(ind,nperm,replace=TRUE)
  }
  SNP_sampling = do.call(cbind,ind_lst)
  TISSUE.eqtl <- paste(strsplit((paste(strsplit(df.to.run$tissue[i],'_-_|\\(|\\)')[[1]],collapse = '_')),'__')[[1]],collapse='_')
  eQTL <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_Analysis_v7_eQTL/',TISSUE.eqtl,'.v7.signif_variant_gene_pairs.txt'),data.table = F)
  # pval.perm.vec <- rep(NA,nperm)
  p.vec <- rep(NA,nperm)
  for (j in 1:nperm) {
    perm_snp = LD$gtex_snpid[SNP_sampling[j,]]
    n = sum(perm_snp %in% eQTL$variant_id)
    N <- nrow(df.tmp)
    p.vec[j] = n/N
    # pval.perm.vec[j] <- binom.test(sum(df.tmp$eQTL),N,p.vec[j])$p.value
  }
  pval.vec[i] <- binom.test(sum(df.tmp$eQTL),N,mean(p.vec))$p.value
  # pval
  # summary(pval.vec)
}


df.to.run$pval <- pval.vec
fwrite(df.to.run,'/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/020819_eQTL_enrichment.matched.csv',quote=F,row.names = F,col.names = T,sep='\t')



# eQTL_enrichment = fread('/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/eQTL_enrichment.1.csv',data.table = F,stringsAsFactors = F)
# head(eQTL_enrichment)
# 
# df.mg <- merge(df.to.run,eQTL_enrichment,by.y=c('Tissue','PhenoNum'),by.x=c('tissue','pheno_num'))
# df.mg$p.value[is.na(df.mg$pval)] <- NA
# wilcox.test(-log10(subset(df.mg,eQTL.enriched==0)$p.value),-log10(subset(df.mg,eQTL.enriched==1)$p.value))
# t.test(-log10(subset(df.mg,eQTL.enriched==0)$p.value),-log10(subset(df.mg,eQTL.enriched==1)$p.value))
# 
# wilcox.test(-log10(subset(df.mg,eQTL.enriched==0)$pval),-log10(subset(df.mg,eQTL.enriched==1)$pval))
# t.test(-log10(subset(df.mg,eQTL.enriched==0)$pval),-log10(subset(df.mg,eQTL.enriched==1)$pval))
# 
# fwrite(df.mg,'/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/012619_eQTL_enrichment.1.csv',quote=F,row.names = F,col.names = T,sep='\t')

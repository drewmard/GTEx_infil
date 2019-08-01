library(data.table)
library(stringr)

# Arguments
# args = commandArgs(trailingOnly=TRUE)
# z <- as.numeric(args[1]) # what pheno to start at?
z <- 1

# init:
nperm=100
tis.old <- ''
num.to.run <- 223
param.df <- matrix(NA,num.to.run,7)
abs <- FALSE

# LD european population patterns - TAKES A WHILE TO READ IN.
print('Reading European LD mapping...')
LD <- fread('/athena/elementolab/scratch/anm2868/GTEx/LD_EUR.tsv',data.table=F,stringsAsFactor=F,sep='\t',header=F)

print('Reading in allele frequencies...')
frq <- fread('/athena/elementolab/scratch/anm2868/GTEx/plink.frq',data.table = F,stringsAsFactors = F)

# gtex-rsid mapping file
print('Reading in number of SNPs in LD...')
gtex_to_rsid <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_rsid_matched.txt',data.table=F,stringsAsFactors=F,header=F)

# number of snps in LD (takes a while)
print('Counting number of SNPs in LD...')
LD$LD_SNP_ct <- str_count(LD[,2],';')+1

# MERGE: rs to gtex names (takes moderate time)
print('Linking rs to gtex SNP ids...')
LD <- merge(LD,gtex_to_rsid,by.x='V1',by.y='V2')

# add in freq 
colnames(LD) <- c('rsid','LD_snps','LD_SNP_ct','gtex_snpid')
print('Merge allele frequencies with LD data...')
LD2 <- merge(LD,frq[,c('SNP','MAF')],by.x='gtex_snpid',by.y='SNP')

# cycle through phenotypes
for (i in z:(z+num.to.run-1)) {
  
  ind_lst <- list()
  
  # What tissue/cell type?
  infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
  tis <- infil_pheno$tissue[i]
  cell <- infil_pheno$cell[i]
  
  if (tis != tis.old) {
    print('Read in cis-eqtl data for that tissue...')
    tis.nospace <- str_replace_all(tis,' ','_')
    TISSUE.eqtl <- paste(strsplit((paste(strsplit(tis.nospace,'_-_|\\(|\\)')[[1]],collapse = '_')),'__')[[1]],collapse='_')
    eQTL <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_Analysis_v7_eQTL/',TISSUE.eqtl,'.v7.signif_variant_gene_pairs.txt'),data.table = F,stringsAsFactors = F)
  }
  
  print(paste0('Reading suggested significant GWAS results for ',i,'...'))
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.sig.txt')
  if (file.exists(f)) {
    df <- fread(f,data.table = F,stringsAsFactors = F)
    df.sub <- subset(df, Pval_Brown < 1e-5)
    if (abs) {
      df.sub <- subset(df, Pval_Brown_Abs < 1e-5)
    }
    
    # merge GWAS to LD results
    df.mg <- merge(df.sub,LD2,by.x='SNP',by.y='gtex_snpid')
    df.mg <- df.mg[,-which(colnames(df.mg)=='LD_snps')] # remove lengthy column
    
    # match snps based on LD and MAF
    ld_diff <- 0
    ind <- c()
    for (j in 1:nrow(df.mg)) {
      while (length(ind) <= 1) {
        ind <- which(LD2$LD_SNP_ct >= (df.mg$LD_SNP_ct[j] - ld_diff) &
                       LD2$LD_SNP_ct <= (df.mg$LD_SNP_ct[j] + ld_diff) &
                       LD2$MAF > (df.mg$MAF[j]-0.01) &
                       LD2$MAF < (df.mg$MAF[j]+0.01))
        ld_diff <- ld_diff + 1
      }
      ind <- ind[ind!=match(df.mg$SNP[j],LD2$gtex_snpid)] # remove identified gwas snp
      ind_lst[[j]] <- sample(ind,nperm,replace=TRUE)
    }
    
    # perform permutation
    SNP_sampling = do.call(cbind,ind_lst)
    p.vec <- rep(NA,nperm)
    eQTL_obs_ct <- sum(df.sub$SNP %in% eQTL$variant_id)
    for (j in 1:nperm) {
      perm_snp = LD2$gtex_snpid[SNP_sampling[j,]]
      n = sum(perm_snp %in% eQTL$variant_id)
      N <- nrow(df.mg)
      p.vec[j] = n/N
    }
    
    # binomial test & save
    p <- binom.test(eQTL_obs_ct,N,mean(p.vec),alternative = "greater")$p.value
    
    # save
    param.df[(i %% num.to.run),] <- c(tis,cell,i,eQTL_obs_ct,N,mean(p.vec),p)
    tis.old <- tis
  } else {
    print(paste0('No significant results: ', i))
  }
}

colnames(param.df) <- c('tissue','cell','i','eQTL_obs_ct','N','mean.p','p')

param.df2 <- data.frame(param.df)
param.df2$p <- (as.numeric(as.character(param.df2$p)))
param.df2$eQTL_obs_ct <- (as.numeric(as.character(param.df2$eQTL_obs_ct)))
param.df2$N <- (as.numeric(as.character(param.df2$N)))
param.df2$mean.p <- (as.numeric(as.character(param.df2$mean.p)))

# two sided test
param.df2$p2 <- NA
for (i in 1:nrow(param.df2)) {
  if (!is.na(param.df2$p[i])) {
    p <- binom.test(param.df2$eQTL_obs_ct[i],param.df2$N[i],param.df2$mean.p[i],alternative = "two.sided")$p.value
    param.df2$p2[i] <- p
  }
}

dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2')
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/eQTL_enrich2_',z,'.txt')
if (abs) {
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/eQTL_enrich2_',z,'.Abs.txt')
}
fwrite(param.df2,f,quote=F,row.names = F,col.names = T,sep='\t',na='NA')

param.df2[order(as.numeric(as.character(param.df2$p))),][1:5,]
param.df2[order(as.numeric(as.character(param.df2$p2))),][1:50,]

# joint across all phenotypes
binom.test(sum(param.df2$eQTL_obs_ct,na.rm = T),sum(param.df2$N,na.rm = T),sum(param.df2$N*param.df2$mean.p,na.rm = T)/sum(param.df2$N,na.rm=T),'greater')$p.value
binom.test(round(mean(param.df2$eQTL_obs_ct,na.rm = T)),round(mean(param.df2$N,na.rm = T)),mean(param.df2$mean.p,na.rm = T),'greater')$p.value
joint_ts_test <- function(x) {binom.test(sum(x$eQTL_obs_ct,na.rm = T),sum(x$N,na.rm = T),sum(x$N*x$mean.p,na.rm = T)/sum(x$N,na.rm=T),'greater')$p.value}
do.call(rbind,lapply(split(param.df2,with(param.df2,tissue),drop=T),joint_ts_test))

# eqtl enriched 2 sided pvalues better than eqtl depleted
param.df2$eQTL.enriched <- as.numeric(param.df2$eQTL_obs_ct/param.df2$N > param.df2$mean.p)
wilcox.test(-log10(subset(param.df2,eQTL.enriched==1)$p2),-log10(subset(param.df2,eQTL.enriched==0)$p2),alternative = 'greater')


# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for combining related p-values using empirical brown's method


# Packages
library(data.table)
library(EmpiricalBrownsMethod)
library(parallel)

# Arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) # what pheno to look at?

# Read in infiltration data
infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
# df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
# df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
# df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
df.pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt',data.table=F,header=T,stringsAsFactors = F)

# initialize
tis <- infil_pheno$tissue[i]
cell <- infil_pheno$cell[i]
df.pheno <- df.pheno[,paste0('pheno',i,'.',1:3)]
df.pheno[df.pheno==-9] <- NA
cellTypes <- c('T cells CD8','CD4_Tcells','Neutrophils','MacrophageSum')
ind <- which(cellTypes %in% cell)
cell2 <- c('CD8Sum','CD4Sum','Neutrophils','MacrophageSum')[ind]

# cibersort rel
print('CIBERSORT (relative)')
df1 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.1.qassoc'),data.table=F,stringsAsFactors = F)
df1 <- df1[,c('CHR','SNP','BP','NMISS','BETA','P')] 
colnames(df1)[(ncol(df1)-1):ncol(df1)] <- c('beta_cibersort_rel','p_lrt_cibersort_rel')

# CIBERSORT (absolute)
print('CIBERSORT (absolute)')
df2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.2.qassoc'),data.table=F,stringsAsFactors = F)
df2 <- df2[,c('SNP','BETA','P')]
colnames(df2)[2:3] <- c('beta_cibersort_abs','p_lrt_cibersort_abs')

df.all <- merge(df1,df2,by='SNP')

# xcell
print('xcell')
df3 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.3.qassoc'),data.table=F,stringsAsFactors = F)
df3 <- df3[,c('SNP','BETA','P')]
colnames(df3)[2:3] <- c('beta_xCell','p_lrt_xCell')

df.all <- merge(df.all,df3,by='SNP')

# arrange for EBM
df.t <- as.data.frame(t(df.all[,c('p_lrt_cibersort_rel',
                               'p_lrt_cibersort_abs',
                               'p_lrt_xCell')]))

# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/issues/1
# Slightly modified version of empiricalBrownsMethod that allows a pre-calculated covariance matrix
# Much more efficient if you have to run empiricalBrownsMethod multiple times with the same data_matrix
empiricalBrownsMethod2 <- function(p_values, extra_info, data_matrix, covar_matrix) {
  if (missing(covar_matrix)) covar_matrix = EmpiricalBrownsMethod:::calculateCovariances(data_matrix)
  return(EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info))
}

# Calculate p-values using all 3 methods:
df.pheno2 <- as.data.frame(t(df.pheno[-which(apply(df.pheno,1,function(x) sum(is.na(x))!=0)),]))
cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno2)
print('Running sped up pre-calculated covariance matrix method...')
p <- unlist(lapply(df.t,function(x) {return(
  empiricalBrownsMethod2(data_matrix = df.pheno,p_values = x,extra_info = F,covar_matrix = cov.matrix))}))
# merge EBM p values
df.all$Pval_Brown <- p

# Calculate p-values using only absolute methods:
df.pheno2 <- as.data.frame(t(df.pheno[-which(apply(df.pheno,1,function(x) sum(is.na(x))!=0)),])[-1,])
cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno2)
print('Running sped up pre-calculated covariance matrix method...')
p <- unlist(lapply(df.t[-1,],function(x) {return(
  empiricalBrownsMethod2(data_matrix = df.pheno,p_values = x,extra_info = F,covar_matrix = cov.matrix))}))
# merge EBM absolute p values
df.all$Pval_Brown_Abs <- p

print('Writing...')  
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.txt')
fwrite(df.all,f,sep='\t',row.names = F,col.names = T,quote = F)

print('Subsetting significant associations...')  
df.sub <- subset(df.all,Pval_Brown < 1e-5 | Pval_Brown_Abs < 1e-5 | p_lrt_cibersort_rel < 1e-5 | p_lrt_cibersort_abs < 1e-5 | p_lrt_xCell < 1e-5 )
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.sig.txt')
fwrite(df.sub,f,sep='\t',row.names = F,col.names = T,quote = F)




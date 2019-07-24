# Packages
library(data.table)
library(EmpiricalBrownsMethod)
library(parallel)

# Arguments
args = commandArgs(trailingOnly=TRUE)
start_ind <- args[1] # what pheno to look at?

# Read in tissues
#tis_list <- fread('/home/anm2868/etc/TISSUES.txt',data.table=F,header=F)[,1]
tis_by_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/df.torun.txt',data.table=F,header=F)

# Initialize
df.results <- data.frame()
df.results2 <- data.frame()

# create save p value directory
dir.create('/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2')

# Iterate through tissues
print('Begin iterating through tissues: ')
rng <- seq(1,nrow(tis_by_pheno),by=3)
for (indicator in start_ind:length(rng)) {
# for (i in c(217)) {  
  i <- rng[indicator]
  tis <- tis_by_pheno[i,1]
  cell <- tis_by_pheno[i,3]
  pheno_num <- tis_by_pheno[i,2]
  # print(paste0('Analysis #',indicator,': ',tis,' (Pheno: ',cell,')'))
  print(paste0('Analysis #',i,': ',tis,' (Pheno: ',cell,')'))
  df.pheno <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name.',tis,'.fam'),data.table=F,header=F)
  df.pheno <- df.pheno[,(5+pheno_num):(5+pheno_num+2)]
  df.pheno <- as.data.frame(t(df.pheno))

  # cibersort rel
  df1 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results/gtex_all.filter.name.',tis,'.LM_out',pheno_num,'.assoc.txt'),data.table=F)
  df1 <- df1[,c(1:10,12)] # remove wald and score test; keep LRT
  colnames(df1)[(ncol(df1)-2):ncol(df1)] <- c('beta_cibersort_rel','se_cibersort_rel','p_lrt_cibersort_rel')

  # CIBERSORT (absolute)
  df2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results/gtex_all.filter.name.',tis,'.LM_out',pheno_num+1,'.assoc.txt'),data.table=F)
  df2 <- df2[,c('rs','beta','se','p_lrt')]
  colnames(df2)[2:4] <- c('beta_cibersort_abs','se_cibersort_abs','p_lrt_cibersort_abs')
  
  df1 <- merge(df1,df2,by='rs')
  
  # xcell
  df2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results/gtex_all.filter.name.',tis,'.LM_out',pheno_num+2,'.assoc.txt'),data.table=F)
  df2 <- df2[,c('rs','beta','se','p_lrt')]
  colnames(df2)[2:4] <- c('beta_xCell','se_xCell','p_lrt_xCell')
  
  df1 <- merge(df1,df2,by='rs')

  # Combine p-values using Brown's method.
  # df.t <- as.data.frame(t(df1[,c('p_lrt_xCell',
  #                            'p_lrt_cibersort_rel',
  #                            'p_lrt_cibersort_abs')]))
  
  df.t <- as.data.frame(t(df1[,c('p_lrt_cibersort_rel',
                              'p_lrt_cibersort_abs',
                              'p_lrt_xCell')]))

  # https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/issues/1
  # Slightly modified version of empiricalBrowsMethod that allows a pre-calculated covariance matrix
  # Much more efficient if you have to run empiricalBrownsMethod multiple times with the same data_matrix
  empiricalBrownsMethod2 <- function(p_values, extra_info, data_matrix, covar_matrix) {
    if (missing(covar_matrix)) covar_matrix = EmpiricalBrownsMethod:::calculateCovariances(data_matrix)
    return(EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info))
  }

  cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno)

  print('Running sped up pre-calculated covariance matrix method...')
  p <- unlist(lapply(df.t,function(x) {return(
    empiricalBrownsMethod2(data_matrix = df.pheno,p_values = x,extra_info = F,covar_matrix = cov.matrix))}))
  df1$Pval_Brown <- p
  
  print('Writing...')  
  fwrite(df1[,c('rs','Pval_Brown')],paste0('/athena/elementolab/scratch/anm2868/GTEx/Brown_P_val_2/',tis,'.',pheno_num,'.txt'),sep='\t',row.names = F,col.names = T,quote = F)
  
}

print('Script completed.')

# Packages
library(data.table)

# Arguments
args = commandArgs(trailingOnly=TRUE)
pheno_num <- args[1] # what pheno to look at?

# Read in tissues
tis_list <- fread('/home/anm2868/etc/TISSUES.txt',data.table=F)[,1]

# Initialize
df.results <- data.frame()
df.results2 <- data.frame()

# Iterate through tissues
print('Begin iterating through tissues: ')
for (i in 1:length(tis_list)) {
  
  tis <- tis_list[i]
  
  print(paste0('Tissue #',i,': ',tis,' (Pheno #',pheno_num,')'))
  
  # xCell
  df1 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results8/gtex_all.filter.name.',tis,'.LM_out',pheno_num,'.assoc.txt'),data.table=F)
  df1 <- df1[,c(1:10,12)] # remove wald and score test; keep LRT
  colnames(df1)[(ncol(df1)-2):ncol(df1)] <- c('beta_xCell','se_xCell','p_lrt_xCell')
  
  # CIBERSORT (relative)
  df2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results7/gtex_all.filter.name.',tis,'.LM_out',pheno_num,'.assoc.txt'),data.table=F)
  df2 <- df2[,c('rs','beta','se','p_lrt')]
  colnames(df2)[2:4] <- c('beta_cibersort_rel','se_cibersort_rel','p_lrt_cibersort_rel')
  
  df1 <- merge(df1,df2,by='rs')
  
  # CIBERSORT (absolute)
  df2 <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/Results11/gtex_all.filter.name.',tis,'.LM_out',pheno_num,'.assoc.txt'),data.table=F)
  df2 <- df2[,c('rs','beta','se','p_lrt')]
  colnames(df2)[2:4] <- c('beta_cibersort_abs','se_cibersort_abs','p_lrt_cibersort_abs')
  
  df1 <- merge(df1,df2,by='rs')
  
  df.tmp <- subset(df1,p_lrt_xCell < 1e-5 & p_lrt_cibersort_rel < 0.05 & p_lrt_cibersort_abs < 0.05)
  # df.tmp <- subset(df1,p_lrt_xCell < 0.05 & p_lrt_cibersort_rel < 1e-5 & p_lrt_cibersort_abs < 0.05)

  if (nrow(df.tmp) > 0) {
    df.tmp$TISSUE <- tis
    if (nrow(df.results) == 0) {
      df.results <- df.tmp
    } else {
      df.results <- rbind(df.results,df.tmp)
    }
  }

  # df.tmp <- subset(df1,p_lrt_xCell < 1e-5 & p_lrt_cibersort_rel < 0.05 & p_lrt_cibersort_abs < 0.05)
  df.tmp <- subset(df1,p_lrt_xCell < 0.05 & p_lrt_cibersort_rel < 1e-5 & p_lrt_cibersort_abs < 0.05)

  if (nrow(df.tmp) > 0) {
    df.tmp$TISSUE <- tis
    if (nrow(df.results2) == 0) {
      df.results2 <- df.tmp
    } else {
      df.results2 <- rbind(df.results2,df.tmp)
    }
  }

}

print('Saving results...')
fwrite(df.results,paste0('/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name.SIGHIT.',pheno_num,'.txt'),sep='\t',row.names = F,col.names = T,quote = F)
fwrite(df.results2,paste0('/athena/elementolab/scratch/anm2868/GTEx/gtex_all.filter.name.SIGHIT.',pheno_num,'.CIBER_rel.txt'),sep='\t',row.names = F,col.names = T,quote = F)


print('Script completed.')

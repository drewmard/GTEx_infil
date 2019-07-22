library('data.table');# library('dplyr')
library('stringr')
df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/xCell_TPM_TBT_mod.txt',data.table = F,stringsAsFactors = F)

df.results <- fread('/athena/elementolab/scratch/anm2868/GTEx/df.torun.txt',data.table=F,stringsAsFactors = F)
df.results <- unique(df.results[,c(1,3)])
colnames(df.results) <- c('tissue','cell')
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)

# Run MLR for age & sex
df.coef <- as.data.frame(matrix(NA,nrow(df.results),14))
for (j in 1:3) {
  for (i in 1:nrow(df.results)) {
    # Subset
    tis.space <- df.results[i,1]
    tis <- str_replace_all(tis.space,'_',' ')
    cell <- df.results[i,2]
    
    # Run cibersort rel
    if (j == 1) {
      df.sub <- subset(df.rel,SMTSD==tis)
      
      df.attr.sub <- subset(df.attr,SMTSD==tis)
      df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]

      df.sub <- merge(df.sub,df.attr.sub,by.x='Input Sample',by.y='SAMPID')

      # Assign NAs a value so that model does not discard data
      df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
      df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4

      # run multiple linear regression
      if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
        # sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY),data=df.sub))$coef
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,3:4] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else if (tis %in% c('Whole Blood')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,3:4] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,5:6] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else {
        # sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX),data=df.sub))$coef
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,3:4] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,5:6] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      }
    }
    
    # Run cibersort absolute
    if (j == 2) {
      df.sub <- subset(df.abs,SMTSD==tis)

      df.attr.sub <- subset(df.attr,SMTSD==tis)
      df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]
      paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
      df.attr.sub$ID <- sapply(strsplit(df.attr.sub$SAMPID,"-"),paste.s)
      df.attr.sub <- df.attr.sub[,-1]
      
      df.sub <- merge(df.sub,df.attr.sub,by='ID')
      
      # Assign NAs a value so that model does not discard data
      df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
      df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4

      # run multiple linear regression
      if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,7:8] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else if (tis %in% c('Whole Blood')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,7:8] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,9:10] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,7:8] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,9:10] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      }
    }
      
    # Run xCell
    if (j == 3) {
      df.sub <- subset(df.xcell,SMTSD==tis)

      cellTypes <- c('T cells CD8','CD4_Tcells','Neutrophils','MacrophageSum','T cells regulatory (Tregs)')
      ind <- which(cellTypes %in% cell)
      # cell <- c('CD8+ T-cells','CD4+ T-cells','Neutrophils','Macrophages','Tregs')[ind]
      cell <- c('CD8Sum','CD4Sum','Neutrophils','MacrophageSum')[ind]
      
      df.attr.sub <- subset(df.attr,SMTSD==tis)
      df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]
      paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
      df.attr.sub$ID <- sapply(strsplit(df.attr.sub$SAMPID,"-"),paste.s)
      df.attr.sub <- df.attr.sub[,-1]
      
      df.sub <- merge(df.sub,df.attr.sub,by.x='ID',by.y='ID')
      
      # Assign NAs a value so that model does not discard data
      df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
      df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4
      
      # run multiple linear regression
      if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,11:12] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else if (tis %in% c('Whole Blood')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,11:12] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,13:14] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,11:12] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,13:14] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      }
    }
      
  }
}
df.coef <- as.data.frame(df.coef)
colnames(df.coef) <- c('tis','pheno',
                       'beta_age_cib.rel','p_age_cib.rel','beta_sex_cib.rel','p_sex_cib.rel',
                       'beta_age_cib.abs','p_age_cib.abs','beta_sex_cib.abs','p_sex_cib.abs',
                       'beta_age_xcell','p_age_xcell','beta_sex_xcell','p_sex_xcell'
                       )

#################################################################

# Merge p-values using Brown's method:
library(EmpiricalBrownsMethod)

# Slightly modified version of empiricalBrowsMethod that allows a pre-calculated covariance matrix
# Much more efficient if you have to run empiricalBrownsMethod multiple times with the same data_matrix
empiricalBrownsMethod2 <- function(p_values, extra_info, data_matrix, covar_matrix) {
  if (missing(covar_matrix)) covar_matrix = EmpiricalBrownsMethod:::calculateCovariances(data_matrix)
  return(EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info))
}

# Initialize
df.coef$p_age_brown <- NA; df.coef$p_sex_brown <- NA
df.coef$age.coef_direc <- NA; df.coef$sex.coef_direc <- NA

for (i in 1:nrow(df.coef)) {
  
  tis.space <- df.results$tissue[i]
  tis <- str_replace_all(tis.space,'_',' ')
  
  # Cibersort rel
  df.sub <- subset(df.rel,SMTSD==tis)
  pheno <- matrix(NA,3,nrow(df.sub))
  cell <- df.results$cell[i]
  pheno[1,] <- df.sub[,cell]
  
  # Cibersort abs
  df.sub <- subset(df.abs,SMTSD==tis)
  pheno[2,] <- df.sub[,cell]
  
  # xCell
  df.sub <- subset(df.xcell,SMTSD==tis)
  ind <- which(cellTypes %in% cell)
  cell <- c('CD8+ T-cells','CD4+ T-cells','Neutrophils','Macrophages','Tregs')[ind]
  pheno[3,] <- df.sub[,cell]
  
  cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(pheno)
  
  # print('Running sped up pre-calculated covariance matrix method...')
  age.p_values <- df.coef[i,c('p_age_cib.rel','p_age_cib.abs','p_age_xcell')]
  ind.min_p <- which.min(age.p_values)
  beta_col <- c('beta_age_cib.rel','beta_age_cib.abs','beta_age_xcell')[ind.min_p]
  beta <- df.coef[i,beta_col]
  beta_direc = NA; if (beta > 0) {beta_direc <- '+'} else {beta_direc <- '-'}
  df.coef$p_age_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno,p_values = age.p_values,extra_info = F,covar_matrix = cov.matrix)
  df.coef$age.coef_direc[i] <- beta_direc
  
  sex.p_values <- df.coef[i,c('p_sex_cib.rel','p_sex_cib.abs','p_sex_xcell')]
  ind.min_p <- which.min(sex.p_values)
  beta_direc = NA; 
  if (length(ind.min_p) > 0) {
    beta_col <- c('beta_sex_cib.rel','beta_sex_cib.abs','beta_sex_xcell')[ind.min_p]
    beta <- df.coef[i,beta_col]
    if (beta > 0) {beta_direc <- '+'} else {beta_direc <- '-'} 
  }
  df.coef$p_sex_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno,p_values = sex.p_values,extra_info = F,covar_matrix = cov.matrix)
  df.coef$sex.coef_direc[i] <- beta_direc
  
}

# Use FDR to correct

df.coef$p_age_brown.fdr <- p.adjust(df.coef$p_age_brown,method='fdr')
df.coef$p_sex_brown.fdr <- p.adjust(df.coef$p_sex_brown,method='fdr')
df.coef$pheno[df.coef$pheno=='CD4Sum'] <- 'CD4+ T cells'
df.coef$pheno[df.coef$pheno=='CD8Sum'] <- 'CD8+ T cells'
df.coef$pheno[df.coef$pheno=='MacrophageSum'] <- 'Macrophages'

df.coef[c(20:22,29:31)]

df.coef.age <- subset(df.coef,p_age_brown.fdr < 0.1)
df.coef.age <- df.coef.age[,c('tis','pheno','age.coef_direc','p_age_brown','p_age_brown.fdr')]
df.coef.age$p_age_brown <- formatC(df.coef.age$p_age_brown, format = "e", digits = 1)
df.coef.age$p_age_brown.fdr <- formatC(df.coef.age$p_age_brown.fdr, format = "e", digits = 1)
colnames(df.coef.age) <- c('Tissue','Cell Type','Effect Direction','P-value','Adj P-value')
library('stargazer')
stargazer(df.coef.age,type='html',summary=F,rownames = F,column.sep.width = "15pt",
          out="/athena/elementolab/scratch/anm2868/GTEx/Age_Assoc_Table.doc")

df.coef.sex <- subset(df.coef,p_sex_brown.fdr < 0.1)
df.coef.sex <- df.coef.sex[,c('tis','pheno','sex.coef_direc','p_sex_brown','p_sex_brown.fdr')]
df.coef.sex$p_sex_brown <- formatC(df.coef.sex$p_sex_brown, format = "e", digits = 1)
df.coef.sex$p_sex_brown.fdr <- formatC(df.coef.sex$p_sex_brown.fdr, format = "e", digits = 1)
colnames(df.coef.sex) <- c('Tissue','Cell Type','Effect Direction','P-value','Adj P-value')
library('stargazer')
stargazer(df.coef.sex,type='html',summary=F,rownames = F,column.sep.width = "15pt",
          out="/athena/elementolab/scratch/anm2868/GTEx/Sex_Assoc_Table.doc")


# Breast vs Cd8 t cell infiltration
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.sub <- subset(df.abs,SMTSD=='Breast - Mammary Tissue')
g <- ggplot(df.sub,aes(x=as.factor(SEX),y=`T cells CD8`,fill=as.factor(SEX))) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.1) + theme_bw() + theme(legend.position='none',panel.grid = element_blank()) +
  labs(x='Sex',y='CD8+ T cells (%)') + scale_x_discrete(labels=c('Male','Female'))
png('/Volumes/SeagateBackupPlusDrive/Elemento/CD8_Sex_Breast.png')
print(g)
dev.off()



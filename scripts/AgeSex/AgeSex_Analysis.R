# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script: analysis of age & sex influences over infiltration patterns

library('data.table')
library('stringr')

# init:
pre_menopause <- FALSE
post_menopause <- FALSE
abs_only <- FALSE

# load infiltration profiles
df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)

# load infiltration phenotypes
infiltration_phenotypes <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)

# load covariate data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)

# cell types data
cellTypes.df <- data.frame( 
  ciber=c('T cells CD8','T cells CD4 naive','CD4_memory','Neutrophils','MacrophageSum',
          'Bcellsum','NK_Sum','DendriticSum','MastSum','Myeloid_Sum',
          'T cells follicular helper','T cells regulatory (Tregs)','T cells gamma delta',
          'Monocytes','Eosinophils','Lymph_Sum'),
  xcell=c('CD8Sum','CD4+ naive T-cells','CD4_memory','Neutrophils','MacrophageSum',
          'Bcellsum','NK cells','DendriticSum','Mast cells','Myeloid_Sum',
          'Th_Sum','Tregs','Tgd cells',
          'Monocytes','Eosinophils','Lymph_Sum'),
  stringsAsFactors = F)


# Run MLR for age & sex
df.coef <- as.data.frame(matrix(NA,nrow(infiltration_phenotypes),14))
for (j in 1:3) {
  for (i in 1:nrow(infiltration_phenotypes)) {
    # Subset
    tis <- infiltration_phenotypes[i,1]
    cell <- infiltration_phenotypes[i,2]
    
    # Run cibersort rel
    if (j == 1) {
      df.sub <- subset(df.rel,SMTSD==tis)
      
      df.attr.sub <- subset(df.attr,SMTSD==tis)
      df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]

      df.sub <- merge(df.sub,df.attr.sub,by.x='Input Sample',by.y='SAMPID')

      # Assign NAs a value so that model does not discard data
      df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
      df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4

      if (pre_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
      } else if (post_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
      }
      
      # run multiple linear regression
      if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,3:4] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else if (tis %in% c('Whole Blood')) {
        sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,3:4] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,5:6] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else {
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

      if (pre_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
      } else if (post_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
      }
      
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

      cell.xcell <- cellTypes.df$xcell[cellTypes.df$ciber==cell]

      df.attr.sub <- subset(df.attr,SMTSD==tis)
      df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]
      paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
      df.attr.sub$ID <- sapply(strsplit(df.attr.sub$SAMPID,"-"),paste.s)
      df.attr.sub <- df.attr.sub[,-1]
      
      df.sub <- merge(df.sub,df.attr.sub,by.x='IID',by.y='ID')
      
      # Assign NAs a value so that model does not discard data
      df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
      df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4
      
      if (pre_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
      } else if (post_menopause) {
        df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
      }
      
      # run multiple linear regression
      if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
        sum.tab <- summary(lm(df.sub[,cell.xcell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,11:12] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else if (tis %in% c('Whole Blood')) {
        sum.tab <- summary(lm(df.sub[,cell.xcell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
        df.coef[i,11:12] <- as.numeric(sum.tab["as.numeric(as.factor(AGE))",c(1,4)])
        df.coef[i,13:14] <- as.numeric(sum.tab["as.factor(SEX)2",c(1,4)])
        df.coef[i,1:2] <- c(tis,cell)
      } else {
        sum.tab <- summary(lm(df.sub[,cell.xcell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + as.factor(SMCENTER),data=df.sub))$coef
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
# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/issues/1
empiricalBrownsMethod2 <- function(p_values, extra_info, data_matrix, covar_matrix) {
  if (missing(covar_matrix)) covar_matrix = EmpiricalBrownsMethod:::calculateCovariances(data_matrix)
  return(EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info))
}

# Initialize
df.coef$p_age_brown <- NA; df.coef$p_sex_brown <- NA
df.coef$age.coef_direc <- NA; df.coef$sex.coef_direc <- NA

for (i in 1:nrow(df.coef)) {
  
  tis <- infiltration_phenotypes$tissue[i]
  cell <- infiltration_phenotypes$cell[i]
  
  # Cibersort rel
  df.sub <- subset(df.rel,SMTSD==tis)
  if (pre_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
  } else if (post_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
  }
  pheno <- matrix(NA,3,nrow(df.sub))
  pheno[1,] <- df.sub[,cell]
  
  # Cibersort abs
  df.sub <- subset(df.abs,SMTSD==tis)
  if (pre_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
  } else if (post_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
  }
  pheno[2,] <- df.sub[,cell]
  
  # xCell
  df.sub <- subset(df.xcell,SMTSD==tis)
  if (pre_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) <= 3)
  } else if (post_menopause) {
    df.sub <- subset(df.sub,as.numeric(as.factor(AGE)) > 3)
  }
  cell.xcell <- cellTypes.df$xcell[cellTypes.df$ciber==cell]
  pheno[3,] <- df.sub[,cell.xcell]
  
  if (cell %in% c('CD4.CD8')) {
    z <- which(apply(pheno,2,function(x) sum(is.na(x))==0))
    pheno <- pheno[,z]
  }
  
  cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(pheno)
  
  # print('Running sped up pre-calculated covariance matrix method...')
  age.p_values <- df.coef[i,c('p_age_cib.rel','p_age_cib.abs','p_age_xcell')]
  ind.min_p <- which.min(age.p_values)
  beta_col <- c('beta_age_cib.rel','beta_age_cib.abs','beta_age_xcell')[ind.min_p]
  beta <- df.coef[i,beta_col]
  beta_direc = NA; if (beta > 0) {beta_direc <- '+'} else {beta_direc <- '-'}
  if (abs_only) {
    df.coef$p_age_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno[-1,],p_values = age.p_values[,-1],extra_info = F,covar_matrix = cov.matrix[,-1])
  } else {
    df.coef$p_age_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno,p_values = age.p_values,extra_info = F,covar_matrix = cov.matrix)
  }
  df.coef$age.coef_direc[i] <- beta_direc
  
  sex.p_values <- df.coef[i,c('p_sex_cib.rel','p_sex_cib.abs','p_sex_xcell')]
  ind.min_p <- which.min(sex.p_values)
  beta_direc = NA; 
  if (length(ind.min_p) > 0) {
    beta_col <- c('beta_sex_cib.rel','beta_sex_cib.abs','beta_sex_xcell')[ind.min_p]
    beta <- df.coef[i,beta_col]
    if (beta > 0) {beta_direc <- '+'} else {beta_direc <- '-'} 
  }
  if (abs_only) {
    df.coef$p_sex_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno[-1,],p_values = sex.p_values[,-1],extra_info = F,covar_matrix = cov.matrix[,-1])
  } else {
    df.coef$p_sex_brown[i] <- empiricalBrownsMethod2(data_matrix = pheno,p_values = sex.p_values,extra_info = F,covar_matrix = cov.matrix)
  }
  df.coef$sex.coef_direc[i] <- beta_direc
  
}

#################################################################
# Use FDR to correct

df.coef$p_age_brown.fdr <- NA
df.coef$p_sex_brown.fdr <- NA
df.coef$p_age_brown.fdr <- p.adjust(df.coef$p_age_brown,method='fdr')
df.coef$p_sex_brown.fdr <- p.adjust(df.coef$p_sex_brown,method='fdr')
# df.coef$pheno[df.coef$pheno=='CD4Sum'] <- 'CD4+ T cells'
# df.coef$pheno[df.coef$pheno=='CD8Sum'] <- 'CD8+ T cells'
# df.coef$pheno[df.coef$pheno=='MacrophageSum'] <- 'Macrophages'

if (pre_menopause) {
  print(1)
  fwrite(df.coef,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_results.pre.txt',sep = '\t',quote=F,row.names = F,col.names = T,na='NA')
} else if (post_menopause) {
  print(2)
  fwrite(df.coef,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_results.post.txt',sep = '\t',quote=F,row.names = F,col.names = T,na='NA')
} else if (abs_only) {
  print(3)
  fwrite(df.coef,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_results.abs.txt',sep = '\t',quote=F,row.names = F,col.names = T,na='NA')
} else {
  print(4)
  fwrite(df.coef,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_results.txt',sep = '\t',quote=F,row.names = F,col.names = T,na='NA')
}

library(data.table)
df.coef <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_results.txt',data.table = F,stringsAsFactors = F)
df.coef$phenotype <- paste(df.coef$tis,df.coef$pheno,sep = '-')

df.coef$p_age_xcell.fdr <- p.adjust(df.coef$p_age_xcell,method='fdr')
df.coef$p_sex_xcell.fdr <- p.adjust(df.coef$p_sex_xcell,method='fdr')
df.coef$p_age_cib.abs.fdr <- p.adjust(df.coef$p_age_cib.abs,method='fdr')
df.coef$p_sex_cib.abs.fdr <- p.adjust(df.coef$p_sex_cib.abs,method='fdr')
df.coef$p_age_cib.rel.fdr <- p.adjust(df.coef$p_age_cib.rel,method='fdr')
df.coef$p_sex_cib.rel.fdr <- p.adjust(df.coef$p_sex_cib.rel,method='fdr')
print(paste0('Age hits: ', nrow(subset(df.coef,p_age_brown.fdr < 0.1))))
print(paste0('Sex hits: ', nrow(subset(df.coef,p_sex_brown.fdr < 0.1))))
print(paste0('Total EBM hits: ', nrow(x <- subset(df.coef,p_age_brown.fdr < 0.1 | p_sex_brown.fdr < 0.1))))
print(paste0('Total xCell hits: ', nrow(subset(df.coef,p_age_xcell.fdr < 0.1 | p_sex_xcell.fdr < 0.1))))
print(paste0('Total CIB-Abs hits: ', nrow(subset(df.coef,p_age_cib.abs.fdr < 0.1 | p_sex_cib.abs.fdr < 0.1))))
print(paste0('Total CIB-Rel hits: ', nrow(subset(df.coef,p_age_cib.rel.fdr < 0.1 | p_sex_cib.rel.fdr < 0.1))))
# print(paste0('Total sig phenotypes: ',nrow(x <- subset(df.coef,p_age_brown.fdr < 0.1 | p_sex_brown.fdr < 0.1))))
print(paste0(length(unique(x$tis)),'/',length(unique(df.coef$tis)),' tissues have > 1 significant phenotype.'))
print('Tissues w/ no significant phenotype:')
unique(df.coef$tis)[which(!(unique(df.coef$tis) %in% unique(x$tis)))]
df.coef[order(df.coef$p_age_brown.fdr),][1:10,]
df.coef[order(df.coef$p_sex_brown.fdr),][1:10,]

subset(df.coef,tis=='Whole Blood')
subset(df.coef,tis=='Thyroid')
subset(df.coef,tis=='Breast - Mammary Tissue')
subset(df.coef,tis %in% c('Nerve - Tibial','Artery - Tibial'))
subset(df.coef,tis %in% c('Nerve - Tibial','Artery - Tibial'))
subset(df.coef,tis %in% c('Artery - Aorta','Artery - Coronary'))


# Save significant associations: 

df.coef.age <- subset(df.coef,p_age_brown.fdr < 0.1)
df.coef.age <- df.coef.age[,c('tis','pheno','age.coef_direc','p_age_brown','p_age_brown.fdr')]
df.coef.age$p_age_brown <- formatC(df.coef.age$p_age_brown, format = "e", digits = 1)
df.coef.age$p_age_brown.fdr <- formatC(df.coef.age$p_age_brown.fdr, format = "e", digits = 1)
colnames(df.coef.age) <- c('Tissue','Cell Type','Effect Direction','P-value','Adj P-value')
library('stargazer')
stargazer(df.coef.age,type='html',summary=F,rownames = F,column.sep.width = "15pt",
          out="/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/Age_Assoc_Table.doc")

df.coef.sex <- subset(df.coef,p_sex_brown.fdr < 0.1)
df.coef.sex <- df.coef.sex[,c('tis','pheno','sex.coef_direc','p_sex_brown','p_sex_brown.fdr')]
df.coef.sex$p_sex_brown <- formatC(df.coef.sex$p_sex_brown, format = "e", digits = 1)
df.coef.sex$p_sex_brown.fdr <- formatC(df.coef.sex$p_sex_brown.fdr, format = "e", digits = 1)
colnames(df.coef.sex) <- c('Tissue','Cell Type','Effect Direction','P-value','Adj P-value')
library('stargazer')
stargazer(df.coef.sex,type='html',summary=F,rownames = F,column.sep.width = "15pt",
          out="/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/Sex_Assoc_Table.doc")


# Breast vs Cd8 t cell infiltration plot
# library(ggplot2)
# workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/'
# df.abs <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
# df.sub <- subset(df.abs,SMTSD=='Breast - Mammary Tissue')
# g <- ggplot(df.sub,aes(x=as.factor(SEX),y=Lymph_Sum,fill=as.factor(SEX))) + geom_boxplot(outlier.shape=NA) +
#   geom_jitter(width=0.1) + theme_bw() + theme(legend.position='none',panel.grid = element_blank()) +
#   labs(x='Sex',y='Lymphocytes (CIBERSORT - Absolute)') + scale_x_discrete(labels=c('Male','Female'))
# png('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/Lymphocyte_Sex_Breast.png')
# print(g)
# dev.off()

# library(ggplot2)
# workdir <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/'
# df.abs <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
# df.rel <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt'),data.table = F,stringsAsFactors = F)
# df.xcell <- fread(paste0(workdir,'GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt'),data.table = F,stringsAsFactors = F)
# 
# tissue <- list('Breast - Mammary Tissue','Thyroid')
# g.plots <- list(); df.mg.melt <- list()
# nplot <- 2
# for (i in 1:nplot) {
#   tis <- tissue[[i]]
#   df.mg <- merge(
#     merge(
#       subset(df.abs,SMTSD==tis)[,c('ID','SEX','Lymph_Sum')],
#       subset(df.rel,SMTSD==tis)[,c('ID','Lymph_Sum')],by='ID'),
#     subset(df.xcell,SMTSD==tis)[,c('IID','Lymph_Sum')],
#     by.x='ID',by.y='IID')
#   colnames(df.mg)[3:5] <- c('CIB-Abs','CIB-Rel','xCell')
#   for (j in 3:5) df.mg[,j] <- rntransform(df.mg[,j])
#   df.mg.melt[[i]] <- (melt(df.mg[,-1],id='SEX'))
#   g.plots[[i]] <- ggplot(df.mg.melt[[i]],aes(x=as.factor(SEX),y=value,fill=as.factor(variable))) + geom_boxplot(outlier.shape=NA) +
#     #geom_jitter(width=0.1) + 
#     theme_bw() + theme(plot.title = element_text(hjust=0.5)) +#+ theme(legend.position='none',panel.grid = element_blank()) +
#     labs(x='Sex',y='Lymphocytes (transformed scores)',title = tis,fill='Method') + scale_x_discrete(labels=c('Male','Female')); #g.plots[[i]]
#   if (i!=nplot) {g.plots[[i]] <- g.plots[[i]] + theme(legend.position='none')}
# }
# library(cowplot)
# plot_grid(g.plots[[1]],g.plots[[2]],ncol=2,rel_widths = c(0.5,0.6))

# df.mg <- merge(
#   merge(
#     subset(df.abs,SMTSD=='Thyroid')[,c('ID','SEX','Lymph_Sum')],
#     subset(df.rel,SMTSD=='Thyroid')[,c('ID','Lymph_Sum')],by='ID'),
#   subset(df.xcell,SMTSD=='Thyroid')[,c('IID','Lymph_Sum')],
#   by.x='ID',by.y='IID')
# colnames(df.mg)[3:5] <- c('CIB-Abs','CIB-Rel','xCell')
# for (i in 3:5) df.mg[,i] <- rntransform(df.mg[,i])
# df.mg.melt <- (melt(df.mg[,-1],id='SEX'))
# g1 <- ggplot(df.mg.melt,aes(x=as.factor(SEX),y=value,fill=as.factor(variable))) + geom_boxplot(outlier.shape=NA) +
#   #geom_jitter(width=0.1) + 
#   theme_bw() +#+ theme(legend.position='none',panel.grid = element_blank()) +
#   labs(x='Sex',y='Lymphocytes (transformed scores)') + scale_x_discrete(labels=c('Male','Female')); g1
# }

# df.sub <- subset(df.abs,SMTSD=='Thyroid')
# g2 <- ggplot(df.sub,aes(x=as.factor(SEX),y=rntransform(Lymph_Sum),fill=as.factor(SEX))) + geom_boxplot(outlier.shape=NA) +
#   geom_jitter(width=0.1) + theme_bw() + theme(legend.position='none',panel.grid = element_blank()) +
#   labs(x='Sex',y='Lymphocytes (CIBERSORT - Absolute)') + scale_x_discrete(labels=c('Male','Female'))
# png('~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/Lymphocyte_Sex_Breast.png')
# print(g2)
# dev.off()

# library(ggplot2)
# df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
# df.sub <- subset(df.abs,SMTSD=='Breast - Mammary Tissue')
# g <- ggplot(df.sub,aes(x=as.factor(SEX),y=`T cells CD8`,fill=as.factor(SEX))) + geom_boxplot(outlier.shape=NA) +
#   geom_jitter(width=0.1) + theme_bw() + theme(legend.position='none',panel.grid = element_blank()) +
#   labs(x='Sex',y='CD8+ T cells score (CIBERSORT - Absolute)') + scale_x_discrete(labels=c('Male','Female'))
# png('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/CD8_Sex_Breast.png')
# print(g)
# dev.off()


# Artery vs Cd4 t cell infiltration plot
# library(ggplot2)
# df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
# df.sub <- subset(df.abs,SMTSD=='Artery - Tibial')
# g <- ggplot(df.sub,aes(x=as.factor(SEX),y=CD4_Tcells,fill=as.factor(AGE))) + geom_boxplot(outlier.shape=NA) +
#   geom_jitter(width=0.1) + theme_bw() + theme(legend.position='none',panel.grid = element_blank()) +
#   labs(x='Age',y='CD4+ T cells score (CIBERSORT - Absolute)') 
# png('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/CD4_Age_ArteryTibial.png')
# print(g)
# dev.off()


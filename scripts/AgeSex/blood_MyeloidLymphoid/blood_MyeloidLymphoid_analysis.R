library('data.table')
library('stringr')

# load infiltration profiles
df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
tis <- 'Whole Blood'
df.sub <- subset(df.rel,SMTSD==tis)
df.sub2 <- subset(df.xcell,SMTSD==tis)

cor.test(df.sub$CD4.CD8,df.sub2$CD4.CD8)
cor.test(df.sub$Myeloid.Lymph,df.sub2$Myeloid.Lymph)


# load covariate data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)

df.attr.sub <- subset(df.attr,SMTSD==tis)
df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]
# paste.s <- function(x) {return(paste(x[1:2],collapse='-'))}
# df.attr.sub$ID <- sapply(strsplit(df.attr.sub$SAMPID,"-"),paste.s)
# df.attr.sub <- df.attr.sub[,-1]
# df.attr.sub <- df.attr.sub[-which(duplicated(df.attr.sub$ID)),]

# df.sub <- merge(df.sub,df.attr.sub,by.x='ID',by.y='ID')
# df.sub2 <- merge(df.sub2,df.attr.sub,by.x='IID',by.y='ID')

df.sub <- merge(df.sub,df.attr.sub,by.x='Input Sample',by.y='SAMPID')
df.sub2 <- merge(df.sub2,df.attr.sub,by.x='SAMP',by.y='SAMPID')

# Assign NAs a value so that model does not discard data
df.sub$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
df.sub$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4
df.sub2$DTHHRDY[is.na(df.sub$DTHHRDY)] <- 5
df.sub2$SMATSSCR[is.na(df.sub$SMATSSCR)] <- 4

# run multiple linear regression
cell <- 'CD4.CD8'
sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
p.cd4.cd8.ciber <- sum.tab[c(
  'as.numeric(as.factor(AGE))'
),4]
cell <- 'Myeloid.Lymph'
sum.tab <- summary(lm(df.sub[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub))$coef
p.myeloid.lymph.ciber <- sum.tab[c(
  'as.numeric(as.factor(AGE))'
),4]
cell <- 'CD4.CD8'
sum.tab <- summary(lm(df.sub2[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub2))$coef
p.cd4.cd8.xcell <- sum.tab[c(
  'as.numeric(as.factor(AGE))'
),4]
cell <- 'Myeloid.Lymph'
sum.tab <- summary(mod <- lm(df.sub2[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub2))$coef
p.myeloid.lymph.xcell <- sum.tab[c(
  'as.numeric(as.factor(AGE))'
),4]
  
####

# for plotting:
cell <- 'Myeloid.Lymph'
mod <- lm(df.sub[,cell]~as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub2,na.action = na.exclude)
summary(lm(residuals(mod)~as.numeric(as.factor(AGE)),data=df.sub))
myeloid.lymph.ciber <- residuals(mod)

mod <- lm(df.sub2[,cell]~as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER),data=df.sub2,na.action = na.exclude)
summary(lm(residuals(mod)~as.numeric(as.factor(AGE)),data=df.sub2))
myeloid.lymph.xcell <- residuals(mod)

fwrite(data.frame(ID=df.sub$ID,myeloid.lymph.xcell,myeloid.lymph.ciber,Age=df.sub$AGE),'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/ratio.df.txt',row.names = F,col.names = T,sep = '\t',na='NA',quote = F)
#############################33

# Merge p-values using Brown's method:
library(EmpiricalBrownsMethod)

# Slightly modified version of empiricalBrowsMethod that allows a pre-calculated covariance matrix
# Much more efficient if you have to run empiricalBrownsMethod multiple times with the same data_matrix
# https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM/issues/1
empiricalBrownsMethod2 <- function(p_values, extra_info, data_matrix, covar_matrix) {
  if (missing(covar_matrix)) covar_matrix = EmpiricalBrownsMethod:::calculateCovariances(data_matrix)
  return(EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info))
}

pheno <- matrix(NA,2,nrow(df.sub))
pheno[1,] <- df.sub[,'Myeloid.Lymph']
pheno[2,] <- df.sub2[,'Myeloid.Lymph']

z <- which(apply(pheno,2,function(x) sum(is.na(x))==0))
pheno <- pheno[,z]

cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(pheno)

age.p_values <- c(0.131,0.009) # hard coded (same as calculated above)

print(empiricalBrownsMethod2(data_matrix = pheno,p_values = age.p_values,extra_info = F,covar_matrix = cov.matrix))

####

pheno <- matrix(NA,2,nrow(df.sub))
pheno[1,] <- df.sub[,'CD4.CD8']
pheno[2,] <- df.sub2[,'CD4.CD8']

z <- which(apply(pheno,2,function(x) sum(is.na(x))==0))
pheno <- pheno[,z]

cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(pheno)

age.p_values <- c(0.005,0.634) # hard coded (same as calculated above)

print(empiricalBrownsMethod2(data_matrix = pheno,p_values = age.p_values,extra_info = F,covar_matrix = cov.matrix))



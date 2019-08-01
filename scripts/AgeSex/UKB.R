library(data.table)
read_data=TRUE

# blood disease individuals
if (read_data) df.disease <- fread('/athena/elementolab/scratch/anm2868/vQTL/UKB/blood_disease.indiv_id.txt',data.table = F,stringsAsFactors = F,header = T)

# UKB sample QC info
if (read_data) df <- fread('/home/kulmsc/athena/ukbiobank/qc/ukb_sample_qc.txt',data.table = F,stringsAsFactors = F)
df2 <- df[,c(4,5,6,19,20,23,24,3,26:45)] # 19=exc heterozyg, 20 = aneuploidy, 23 = exc rel,24 = brits,3= genotyping array
colnames(df2)[1:8] <- c(
  'batch',
  'plate',
  'well',
  'het.missing.outliers',
  'putative.sex.chromosome.aneuploidy',
  'excess.relatives',
  'in.white.British.ancestry.subset',
  'genotyping.array')
colnames(df2)[9:ncol(df2)] <- paste0('PC',1:20)

# phenotype file (covariates)
if (read_data) pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
pheno2 <- pheno[,c('eid','22001-0.0','21022-0.0','22000-0.0','22007-0.0','22008-0.0','22009-0.1',
                   '21001-0.0')]
colnames(pheno2)[2:ncol(pheno2)] <- c('sex','age','batch','plate','well','p.c.1',
                                      'bmi')
# Phenotype file 2 (covariates)
if (read_data) pheno.new <- fread('/home/kulmsc/athena/ukbiobank/setup_morePhenos/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)
# pheno.new2 <- pheno.new[,c('eid',paste0('X',c('3700.0.0','2724.0.0','1558.0.0','1239.0.0','1249.0.0')))]
pheno.new2 <- pheno.new[,c('eid','3700-0.0','2724-0.0','1558-0.0','1239-0.0','1249-0.0')]
colnames(pheno.new2)[2:ncol(pheno.new2)] <- c('time.since.period','menopause','alcohol.freq','current.smoking','past.smoking')
pheno.new2 <- subset(pheno.new2, !(eid %in% df.disease$eid))

# phenotypes
pheno3 <- pheno[,c('eid','30120-0.0','30130-0.0','30140-0.0','30180-0.0','30190-0.0','30200-0.0')]
pheno3 <- subset(pheno3, !(eid %in% df.disease$eid))
PHENOTYPE_NAMES <- c('lymphocyte.count','monocyte.count','neutrophil.count','lymphocyte.percentage','monocyte.percentage','neutrophil.percentage')
colnames(pheno3)[2:ncol(pheno3)] <- PHENOTYPE_NAMES

# Read in fam & merge w/ covariate & phenotype data
if (read_data) fam <- fread('/home/kulmsc/athena/ukbiobank/calls/ukbb.1.fam',data.table = F,stringsAsFactors = F)
fam2 <- cbind(fam,df2)

# Neale subset:
if (read_data) Neale_subset <- read.table('/athena/elementolab/scratch/anm2868/vQTL/UKB/Neale_GWAS/samples.both_sexes.tsv.bgz',header=T,stringsAsFactors = F)
Neale_subset$In <- 1

# Merge diff files together
phenotypeDataFile1 <- merge(pheno2,pheno3,by='eid')
phenotypeDataFile1 <- merge(phenotypeDataFile1,pheno.new2,by='eid')
fam2 <- merge(fam2,phenotypeDataFile1,by.x='V1',by.y='eid',all.x=TRUE)
fam2 <- merge(fam2,Neale_subset,all.x=TRUE,by.x=c('plate.x','well.x'),by.y=c('plate_name','well'))
colnames(fam2)[colnames(fam2) %in% c('V1','V2')] <- c('FID','IID')

# assign NA phenotypes to individuals not in Neale subset or don't pass sample QC
i.remove <- which(fam2$excess.relatives==1 |
                    fam2$putative.sex.chromosome.aneuploidy==1 |
                    fam2$het.missing.outliers==1)
fam2$QC_In <- 1; fam2$QC_In[i.remove] <- 0
fam2$In[which(is.na(fam2$In))] <- 0
for (i in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[i]
  fam2[,paste0(phenoName,'.na')] <- fam2[,phenoName]
  fam2[,paste0(phenoName,'.na')][which(fam2$In==0)] <- NA
  fam2[,paste0(phenoName,'.na')][which(fam2$QC_In==0)] <- NA
}

# remove individuals w/ blood disease
fam2 <- subset(fam2, !(IID %in% df.disease$eid))

fam2.80 <- fam2

# impute age, sex
fam2.80$sex[which(is.na(fam2.80$sex))] <- median(fam2.80$sex,na.rm=T)
# square sex
fam2.80$age2 <- fam2.80$age^2

# impute bmi
tmp <- fam2.80$bmi
tmp[which(abs(scale(tmp)) > 5)] <- NA
val <- median(tmp,na.rm=T)
fam2.80$bmi2 <- fam2.80$bmi;
fam2.80$bmi2.dummy <- as.numeric(is.na(fam2.80$bmi))
i <- which(is.na(fam2.80$bmi2))
j <- which(abs(scale(fam2.80$bmi)) > 5)
fam2.80$bmi2[i] <- val
fam2.80$bmi2.with_outliers <- fam2.80$bmi2
fam2.80$bmi2[j] <- NA

# other variables:
fam2.80$Smoking <- NA
fam2.80$Smoking[which( (fam2.80$current.smoking==1 | fam2.80$current.smoking==2) | (fam2.80$past.smoking==1 | fam2.80$past.smoking==2) )] <- 1
fam2.80$Smoking[which( fam2.80$current.smoking==0 & (fam2.80$past.smoking==3 | fam2.80$past.smoking==4) )] <- 0
fam2.80$Smoking.dummy <- as.numeric(is.na(fam2.80$Smoking))
fam2.80$Smoking[which(is.na(fam2.80$Smoking))] <- median(fam2.80$Smoking,na.rm=T)

tmp <- fam2.80$time.since.period
tmp[which(abs(scale(tmp)) > 5)] <- NA
val <- median(tmp,na.rm=T)
fam2.80$time.since.period2 <- fam2.80$time.since.period;
fam2.80$time.since.period2[which(fam2.80$time.since.period2 %in% c(-1,-3))] <- NA
fam2.80$time.since.period2.dummy <- as.numeric(is.na(fam2.80$time.since.period2))
i <- which(is.na(fam2.80$time.since.period2))
j <- which(abs(scale(fam2.80$time.since.period2)) > 5)
fam2.80$time.since.period2[i] <- val
fam2.80$time.since.period2.with_outliers <- fam2.80$time.since.period2
fam2.80$time.since.period2[j] <- NA

fam2.80$alcohol.freq2 <- fam2.80$alcohol.freq;
fam2.80$alcohol.freq2[which(fam2.80$alcohol.freq %in% c(-3))] <- NA
fam2.80$alcohol.freq2.dummy <- as.numeric(is.na(fam2.80$alcohol.freq2))
fam2.80$alcohol.freq2[which(is.na(fam2.80$alcohol.freq2))] <- median(fam2.80$alcohol.freq2,na.rm=T)

fam2.80$menopause2 <- fam2.80$menopause
fam2.80$menopause2[which(fam2.80$menopause==-3)] <- NA
fam2.80$menopause2[which(fam2.80$menopause==3)] <- NA
fam2.80$menopause2[which(is.na(fam2.80$menopause2))] <- 3
fam2.80$menopause2 <- as.factor(fam2.80$menopause2)

# remove phenotype outliers
for (k in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[k]
  i.outlier <- which(abs(scale(fam2.80[,phenoName])) > 5)
  fam2.80[,paste0(phenoName,'.na')][i.outlier] <- NA
}

# dataframe w/ outliers
fam3.80 <- fam2.80
fam3.80$bmi2 <- fam3.80$bmi2.with_outliers
fam3.80$time.since.period2 <- fam3.80$time.since.period2.with_outliers

# remove phenotype outliers
source('/home/anm2868/scripts/Useful_scripts/rntransform.R')
for (k in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[k]
  suffix<-'.rint'; fam2.80[,paste0(phenoName,'.na',suffix)] <- rntransform(fam2.80[,paste0(phenoName,'.na')])
  suffix<-'.log'; fam2.80[,paste0(phenoName,'.na',suffix)] <- log10(fam2.80[,paste0(phenoName,'.na')]+1)
  
  suffix<-'.rint'; fam3.80[,paste0(phenoName,'.na',suffix)] <- rntransform(fam3.80[,paste0(phenoName,'.na')])
  suffix<-'.log'; fam3.80[,paste0(phenoName,'.na',suffix)] <- log10(fam3.80[,paste0(phenoName,'.na')]+1)
}

df.coef.age <- data.frame(pheno=PHENOTYPE_NAMES,var='age',b=NA,p=NA)
df.coef.sex <- data.frame(pheno=PHENOTYPE_NAMES,var='sex',b=NA,p=NA)
for (k in 1:length(PHENOTYPE_NAMES)) {
  phenoName <- PHENOTYPE_NAMES[k]
  print(paste0('phenotype: ',phenoName))
  
  suffix <- ''
  # male
  mod.formula.1 <- formula(paste(paste0(phenoName,'.na',suffix),' ~ sex+age*sex+genotyping.array+
                               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                               PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+
                               Smoking+Smoking.dummy+
                               alcohol.freq2+alcohol.freq2.dummy+bmi2.dummy+bmi2+bmi2*age'))
  # fitting models, not including individuals that are outliers
  mod1 <- lm(mod.formula.1,
             data=fam2.80,na.action=na.exclude)
  df.coef.age[k,c('b','p')] <- summary(mod1)$coef['age',c('Estimate','Pr(>|t|)')]
  df.coef.sex[k,c('b','p')] <- summary(mod1)$coef['sex',c('Estimate','Pr(>|t|)')]
}

df.coef <- rbind(df.coef.age,df.coef.sex)
fwrite(df.coef,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/AgeSex/AgeSex_UKB_results.txt',sep = '\t',quote=F,row.names = F,col.names = T,na='NA')



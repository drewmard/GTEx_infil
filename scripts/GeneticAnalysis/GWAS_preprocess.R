# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script to prepare data for GWAS

library('data.table')
library('stringr')
source('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/rntransform.R')

############################################################
############################################################
############################################################

# load data
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F)
df.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table = F,stringsAsFactors = F)
tissue_x_celltype <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
fam <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/bin/gtex_all.filter.name.fam',data.table = F,stringsAsFactors = F)

# function that takes as input a tissue type, and outputs the filename to the tissue of interest
tissue_file_name <- function(x.tis) {paste0('/athena/elementolab/scratch/anm2868/GTEx/GENO_PCA/gtex_all.filter.name.',x.tis,'.eigenvec')}

# initialize
colnames(df.xcell)[colnames(df.xcell)=='SAMP'] <- 'Input Sample'
colnames(df.xcell)[colnames(df.xcell)=='IID'] <- 'ID'
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

f <- '/athena/elementolab/scratch/anm2868/GTEx/GENO_PCA/gtex_all.filter.name.eigenvec'
pca <- fread(f,data.table=F,stringsAsFactors = F)
pca <- pca[,-1]; colnames(pca) <- c('ID',paste0('PC',1:min(10,ncol(pca)-1)))

# cycle through infiltration phenotypes:
for (i in 1:nrow(tissue_x_celltype)) {
  tis <- tissue_x_celltype[i,1]
  cell <- tissue_x_celltype[i,2]
  
  # Read in tissue-specific PCA matrix: use first 10 PCs
  # f.tis <- str_replace_all(tissue_file_name(tis),' ','_')
  # pca <- fread(f.tis,data.table=F,stringsAsFactors = F)
  # pca <- pca[,-1]; colnames(pca) <- c('ID',paste0('PC',1:min(10,ncol(pca)-1)))
  
  for (k in 1:3) {
    
    # Indicate which deconvolution outputs to use with k variable
    if (k==1) {
      df <- df.rel
    } else if (k==2) {
      df <- df.abs
    } else if (k==3) {
      df <- df.xcell
      cell <- cellTypes.df$xcell[cellTypes.df$ciber==cell]
    }
    
    # Merging data frames
    df.sub <- subset(df,SMTSD==tis)
    df.sub.pca <- merge(df.sub,pca,by='ID')
    df.attr.sub <- subset(df.attr,SMTSD==tis)
    df.attr.sub <- df.attr.sub[,c('SAMPID','SMATSSCR','SMCENTER')]
    df.sub.pca <- merge(df.sub.pca,df.attr.sub,by.x='Input Sample',by.y='SAMPID')
    
    # Assign NAs a value so that model does not discard data
    df.sub.pca$DTHHRDY[is.na(df.sub.pca$DTHHRDY)] <- 5
    df.sub.pca$SMATSSCR[is.na(df.sub.pca$SMATSSCR)] <- 4
    
    # Run multiple regression to regress out covariates
    if (tis %in% c('Ovary','Uterus','Vagina','Testis','Prostate')) {
      mod <- lm(df.sub.pca[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.numeric(SMATSSCR) + as.factor(SMCENTER) + PC1 + PC2 + PC3,data=df.sub.pca,na.action=na.exclude)
    } else if (tis %in% c('Whole Blood')) {
      mod <- lm(df.sub.pca[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.factor(SMCENTER) + PC1 + PC2 + PC3,data=df.sub.pca,na.action=na.exclude)
    } else {
      mod <- lm(df.sub.pca[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + as.factor(SMCENTER) + PC1 + PC2 + PC3,data=df.sub.pca,na.action=na.exclude)
        # mod <- lm(df.sub.pca[,cell]~as.numeric(as.factor(AGE)) + as.factor(DTHHRDY) + as.factor(SEX) + as.numeric(SMATSSCR) + PC1 + PC2 + PC3,data=df.sub.pca)
    }
    
    # Saving residuals
    resid <- residuals(mod)
    resid.rint <- rntransform(resid) # RINT residuals

    # Save to fam file
    pheno <- data.frame(ID=df.sub.pca$ID,PHENO=resid.rint,stringsAsFactors = F)
    col_name <- paste0('pheno',i,'.',k)
    colnames(pheno)[2] <- col_name;
    fam <- merge(fam,pheno,by.x='V1',by.y='ID',all.x=TRUE)

  }
}
# the missing phenotype
fam <- fam[,-(3:6)]
colnames(fam)[1:2] <- c('FID','IID')
fam[is.na(fam)] <- -9
dir.create('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis',showWarnings = F)
fwrite(fam,'/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt',sep = '\t',col.names = T,row.names = F,quote=F)


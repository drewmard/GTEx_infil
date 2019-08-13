# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to pull all ieGenes from ieQTLs and output for a GeneMania network

library(data.table)
library(stringr)

tis.old <- ''
eQTL2.save <- NA
FIRST <- TRUE
abs <- FALSE

infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes2.txt',data.table=F,stringsAsFactors = F)
infil_pheno <- subset(infil_pheno,!(cell %in% c('CD4.CD8','Myeloid.Lymph')))

# infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)
# infil_pheno <- infil_pheno[1:223,]
# infil_pheno <- subset(infil_pheno,!((tissue=='Heart - Atrial Appendage' & cell=='NK_Sum') | (tissue=='Colon - Transverse' & cell=='NK_Sum') | (tissue=='Cells - EBV-transformed lymphocytes')))

for (i in 1:nrow(infil_pheno)) {
  
  ind_lst <- list()
  
  # What tissue/cell type?
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
    df.sub$tissue <- tis
    df.sub$cell <- cell
    
    eQTL2 <- merge(df.sub,eQTL[,c('variant_id','gene_id')],by.x='SNP',by.y='variant_id')
    
    if (FIRST) {
      eQTL2.save <- eQTL2
      FIRST <- FALSE
    } else {
      eQTL2.save <- rbind(eQTL2.save,eQTL2)
    }
    
    tis.old <- tis
  } else {
    print(paste0('No significant results: ', i))
  }
}

f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.txt')
if (abs) {
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.Abs.txt')
}
fwrite(eQTL2.save,f,sep='\t',quote=F,row.names = F,col.names = T)

f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.ieGene.txt')
if (abs) {
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.ieGene.Abs.txt')
}
x <- do.call(rbind,lapply(strsplit(unique(eQTL2.save$gene_id),'\\.'),function(x) x[1]))
fwrite(as.data.frame(x),f,sep='\t',quote=F,row.names = F,col.names = F)



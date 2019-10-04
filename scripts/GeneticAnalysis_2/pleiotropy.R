# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to assess whether genetic effects exist across cell types and tissues

library(data.table)
library(stringr)

iqtl.save <- NA
FIRST <- TRUE
abs <- FALSE

infil_pheno <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/infiltration_phenotypes.txt',data.table=F,stringsAsFactors = F)

for (i in 1:nrow(infil_pheno)) {
  
  ind_lst <- list()
  
  # What tissue/cell type?
  tis <- infil_pheno$tissue[i]
  cell <- infil_pheno$cell[i]
  
  print(paste0('Reading suggested significant GWAS results for ',i,'...'))
  f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/GWAS/GTEx.pheno',i,'.ALL_EBM.sig.txt')
  if (file.exists(f)) {
    df <- fread(f,data.table = F,stringsAsFactors = F)
    df.sub <- subset(df, Pval_Brown < 1e-5)
    # if (abs) {
    #   df.sub <- subset(df, Pval_Brown_Abs < 1e-5)
    # }
    df <- df.sub
    df.sub$tissue <- tis
    df.sub$cell <- cell
    
    if (FIRST) {
      iqtl.save <- df.sub
      FIRST <- FALSE
    } else {
      iqtl.save <- rbind(iqtl.save,df.sub)
    }
    
  } else {
    print(paste0('No significant results: ', i))
  }
}

# iqtl summary results for supplementary file
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_iQTL.txt')
fwrite(iqtl.save,f,sep='\t',quote=F,row.names = F,col.names = T)

# iQTLs in multiple tissues & save supplementary file
iqtl.save <- fread('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_iQTL.txt',data.table = F,stringsAsFactors = F)
x <- iqtl.save[,c('SNP','tissue')]
x <- x[-which(duplicated(x)),]
y <- sort(table(x$SNP),decreasing = T)
y <- names(y[y>1])
table(table(x$SNP)); table(table(x$SNP))[1]/length(unique(x$SNP))
tab <- table(x$SNP); x <- subset(iqtl.save,SNP%in%names(tab[tab>1])); x <- x[order(x$SNP),]
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_iQTL_MultiTissue.txt')
fwrite(x,f,sep='\t',quote=F,row.names = F,col.names = T,na='NA')

# iQTLs in multiple cell types
x <- iqtl.save[,c('SNP','cell')]
x <- x[-which(duplicated(x)),]
y <- sort(table(x$SNP),decreasing = T)
y <- names(y[y>1])
table(table(x$SNP)); table(table(x$SNP))[1]/(length(unique(x$SNP)))
sort(table(x$SNP),decreasing = T)[1:5]

iqtl.save <-  subset(iqtl.save,Pval_Brown < 5e-8)
aggregate(iqtl.save$Pval_Brown,list(iqtl.save$tissue,iqtl.save$cell),min)
table(table(iqtl.save$SNP))

f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.txt')
eQTL2.save <- fread(f,data.table = F,stringsAsFactors = F)
x <- eQTL2.save[,c('gene_id','tissue')]
x <- x[-which(duplicated(x)),]
y <- sort(table(x$gene_id),decreasing = T)
table(table(x$gene_id))
table(table(x$gene_id))[1]/length(unique(x$gene_id))

# iegenes in multiple tissues & save supplementary file
y <- names(y[y>1])
z <- subset(eQTL2.save[,c('SNP','gene_id','tissue','cell','Pval_Brown')],gene_id %in% y)
z <- z[-which(duplicated(z[,c('gene_id','tissue')])),]
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis2/GWAS_eQTL-only.ieGene.multi.txt')
z2 <- do.call(rbind,lapply(strsplit(unique(z$gene_id),'\\.'),function(x) x[1]))
fwrite(as.data.frame(z2),f,sep='\t',quote=F,row.names = F,col.names = F)

# iegenes in multiple cell types
x <- eQTL2.save[,c('gene_id','cell')]
x <- x[-which(duplicated(x)),]
y <- sort(table(x$gene_id),decreasing = T)
table(table(x$gene_id))
table(table(x$gene_id))[1]/length(unique(x$gene_id))

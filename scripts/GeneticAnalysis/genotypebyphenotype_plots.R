# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# script to draw genotype by phenotype plots

# Arguments
args = commandArgs(trailingOnly=TRUE)
snp <- as.numeric(args[1]) # what pheno to look at?
pheno <- as.numeric(args[2])
# snp <- '1_751756_T_C_b37'

f <- '/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt'
df.pheno <- fread(f,data.table = F,stringsAsFactors = F)
f <- paste0('/athena/elementolab/scratch/anm2868/GTEx/GTEx_infil/output/GeneticAnalysis/Genotype_Subset/gtex_all.',snp,'.raw')
df.geno <- fread(f,stringsAsFactors = F,data.table = F)

df <- merge(df.pheno[,c(2,pheno)],df.geno[,c(2,7)],by='IID')
df <- subset(df,!is.na(df[,ncol(df)]) | (df[,pheno] != -9))
g1 <- ggplot(df,aes(x=as.factor(df[,ncol(df)]),y=df[,2],fill=as.factor(df[,ncol(df)]))) + geom_boxplot(na.rm = T,outlier.shape = NA) + geom_jitter(width=0.1,alpha=0.1) +
  labs(x='rsid (Gene name)',y='Immune cell type (residuals)',title='Tissue type') +
  theme(legend.position="none") +
  # scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))

# save:
scale <- 3
tiff(paste0('/Volumes/SeagateBackupPlusDrive/Elemento/120418/','Geno_x_Pheno','.png'),width=800*scale,height=400*scale,res=100*scale,units="px")
plot_grid(g1,ncol=1)
dev.off()


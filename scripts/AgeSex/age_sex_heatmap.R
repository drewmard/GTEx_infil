library(data.table)
df.coef <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/AgeSex/AgeSex_results.txt',data.table = F,stringsAsFactors = F)
# subset(df.coef,tis %in% c('Nerve - Tibial','Artery - Tibial'))[,c('tis','pheno','p_age_brown','p_age_brown.fdr')]

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
celltype.label.vec <- c('CD8+ T cells','CD4+ naive T cells','CD4+ memory T cells','Neutrophils','Macrophages',
                        'B cells','NK cells','Dendritic cells','Mast cells','Myeloid cells',
                        'Helper T cells','Tregs','Gamma delta T cells',
                        'Monocytes','Eosinophils','Lymphocytes')
df.coef$cell <- celltype.label.vec[match(df.coef$pheno,cellTypes.df$ciber)]

df.coef$Age_Sig <- 0
df.coef$Age_Sig[df.coef$p_age_brown.fdr < 0.1] <- 1
df.coef$Age_Sig[df.coef$p_age_brown.fdr < 0.01] <- 2
df.coef$Age_Sig[df.coef$p_age_brown.fdr < 0.001] <- 3
df.coef$Age_Sig[df.coef$p_age_brown.fdr < 0.00001] <- 4
df.coef$age.coef_direc[df.coef$Age_Sig==0] <- ''

library(ggplot2)
g1 <- ggplot(df.coef,aes(tis,cell)) + 
  geom_tile(aes(fill=as.factor(Age_Sig)),col='white') +
  geom_text(aes(label=age.coef_direc)) +
  scale_fill_manual(values=c('grey','lightgreen','steelblue3','orange','orangered'),
                    labels=c('> 0.1','< 0.1','< 0.01','< 0.001','< 0.00001')) +
  theme_minimal() + theme(panel.grid=element_blank(),
                          # legend.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          # axis.text.x = element_blank()) +
                          axis.text.x = element_text(angle = 45, hjust = 1,size = rel(1))) +
  labs(x='Tissue',y='Cell type',fill='FDR')

df.coef$Sex_Sig <- 0
df.coef$Sex_Sig[df.coef$p_sex_brown.fdr < 0.1] <- 1
df.coef$Sex_Sig[df.coef$p_sex_brown.fdr < 0.01] <- 2
df.coef$Sex_Sig[df.coef$p_sex_brown.fdr < 0.001] <- 3
df.coef$Sex_Sig[df.coef$p_sex_brown.fdr < 0.00001] <- 4
df.coef$sex.coef_direc[df.coef$Sex_Sig==0] <- ''

library(ggplot2)
g2 <- ggplot(df.coef,aes(tis,cell)) + 
  geom_tile(aes(fill=as.factor(Sex_Sig)),col='white') +
  geom_text(aes(label=sex.coef_direc)) +
  scale_fill_manual(values=c('grey','lightgreen','steelblue3','orange','orangered'),
                    labels=c('> 0.1','< 0.1','< 0.01','< 0.001','< 0.00001')) +
  theme_minimal() + theme(panel.grid=element_blank(),
                          axis.text.x = element_text(angle = 45, hjust = 1,size = rel(1))) +
  labs(x='Tissue',y='Cell type',fill='FDR')

library(cowplot)
plot_grid(g2,g1,ncol=2,rel_widths = c(1.275,1))


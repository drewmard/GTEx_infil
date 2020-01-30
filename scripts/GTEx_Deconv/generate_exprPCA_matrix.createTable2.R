f<-'/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/Expr_PCA/perc.var.exp.df.save.txt'
perc.var.exp.df.save <- fread(f,data.table = F,stringsAsFactors = F)
f <- '~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/PCA_infil_cor_full.txt'
PCA_infil_cor_full <- fread(f,data.table = F,stringsAsFactors = F)

perc.var.exp.df.save$PCnum <- substring((perc.var.exp.df.save$PC),3)

PCA_infil_cor_full$ID <- 1:nrow(PCA_infil_cor_full)
df.mg <- merge(PCA_infil_cor_full,perc.var.exp.df.save[,c('Tissue','Percent_Variance_Explained','PCnum')],by.x=c('tissue','PCnum'),by.y=c('Tissue','PCnum'))
df.mg <- df.mg[order(df.mg$ID),]
df.mg <- df.mg[,-which(colnames(df.mg) %in% c('ID'))]

fwrite(df.mg,'~/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GTEx_Deconv/PCA_infil_cor_full.with_perc_var_exp.txt',,quote = F,row.names = F,sep='\t',col.names = T,na='NA')
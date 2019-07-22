i=5 # choose i anywhere between 1 and 53 to select a different tissue
df.attr <- fread('/athena/elementolab/scratch/anm2868/GTEx/COVAR/GTEx_v7_Annotations_SampleAttributesDS.txt',data.table = F,stringsAsFactors = F)
tis.uniq = unique(df.attr$SMTSD)
tis = tis.uniq[i]
print(tis) # What tissue
xcell.sub <- fread(paste0('/athena/elementolab/scratch/anm2868/GTEx/infil_output/xcell_sub_out/','XCell_',tis,'.txt'),data.table=F,stringsAsFactors=F)
df.ciber <- fread('/athena/elementolab/scratch/anm2868/GTEx/infil_output/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table=F,stringsAsFactors=F)
ciber.sub <- subset(df.ciber,SMTSD==tis)

# individuals where no cell type estimated but subtypes are? weird.
subset(xcell.sub,Macrophages==0)[,c('Macrophages M1','Macrophages M2')]

xcell.sub$MacrophageSum <- apply(xcell.sub[,32:34],1,sum)
xcell.sub$CD4Sum <- apply(xcell.sub[,7:10],1,sum)
xcell.sub$CD8Sum <- apply(xcell.sub[,11:14],1,sum)

cor(xcell.sub$Macrophages,ciber.sub$MacrophageSum)
cor(xcell.sub$MacrophageSum,ciber.sub$MacrophageSum)
cor(xcell.sub$`CD4+ T-cells`,ciber.sub$CD4_Tcells)
cor(xcell.sub$CD4Sum,ciber.sub$CD4_Tcells)
cor(xcell.sub$`CD8+ T-cells`,ciber.sub$`T cells CD8`)
cor(xcell.sub$CD8Sum,ciber.sub$`T cells CD8`)

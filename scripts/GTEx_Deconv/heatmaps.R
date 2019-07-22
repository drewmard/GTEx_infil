# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for heatmaps, supp fig 1

library(data.table)
library(ggplot2)

# required files:
df.xcell <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/XCell.all_tissues.txt',data.table=F,stringsAsFactors = F)
df.ciber.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)
df.ciber.abs <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/infiltration_profiles/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt',data.table = F,stringsAsFactors = F)

# RELATIVE:

# subset only tissues w N >= 70 samples
tis.count <- data.frame(table(df.ciber.rel$SMTSD))
df.ciber.rel <- subset(df.ciber.rel,SMTSD %in% subset(tis.count,Freq >= 70)[,1])

# pre-process infil file for plotting:
x <- aggregate(df.ciber.rel[,c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')],list(Tissue=df.ciber.rel$SMTSD),mean)
colnames(x) <- c('Sample Type','CD8+ T cells','CD4+ T cells','Macrophages','Neutrophils')
x1 <- melt(x)

# plot heatmap:
g1 <- ggplot(x1,aes(reorder(`Sample Type`,value,mean,order=TRUE),variable,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       name="Mean\nAbundance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 70,
                                   size = rel(0.85), hjust = 1),
        axis.text.y = element_text(size = rel(0.85)),
        plot.title = element_text(hjust=0.5)) +
  coord_fixed() +
  labs(x='Sample Type',y='Score',title='CIBERSORT - Relative') #+
  # theme(axis.title.x=element_blank(),
  #   axis.text.x=element_blank(),
  #   axis.ticks.x=element_blank())



# ABSOLUTE:

# subset only tissues w N >= 70 samples
tis.count <- data.frame(table(df.ciber.abs$SMTSD))
df.ciber.abs <- subset(df.ciber.abs,SMTSD %in% subset(tis.count,Freq >= 70)[,1])

# pre-process infil file for plotting:
x <- aggregate(df.ciber.abs[,c('T cells CD8','CD4_Tcells','MacrophageSum','Neutrophils')],list(Tissue=df.ciber.abs$SMTSD),mean)
colnames(x) <- c('Sample Type','CD8+ T cells','CD4+ T cells','Macrophages','Neutrophils')
x2 <- melt(x)

# plot heatmap:
g2 <- ggplot(x2,aes(reorder(`Sample Type`,value,mean,order=TRUE),variable,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       name="Mean\nAbundance") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 70,
                                   size = rel(0.85), hjust = 1),
        axis.text.y = element_text(size = rel(0.85)),
        plot.title = element_text(hjust=0.5)) +
  coord_fixed() +
  labs(x='Sample Type',y='Score',title='CIBERSORT - Absolute') #+
# theme(axis.title.x=element_blank(),
#   axis.text.x=element_blank(),
#   axis.ticks.x=element_blank())


# xCell:

# subset only tissues w N >= 70 samples
tis.count <- data.frame(table(df.xcell$SMTSD))
df.xcell <- subset(df.xcell,SMTSD %in% subset(tis.count,Freq >= 70)[,1])
df <- df.xcell

x3 <- aggregate(df[,c('CD8Sum','CD4Sum','MacrophageSum','Neutrophils')],list(Tissue=df$SMTSD),mean)
colnames(x3) <- c('Sample Type','CD8+ T cells','CD4+ T cells','Macrophages','Neutrophils')
x3 <- melt(x3)
g3 <- ggplot(x3,aes(reorder(`Sample Type`,value,mean,order=TRUE),variable,fill=value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       name="Mean\nAbundance") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 70, #vjust = 1,
                                   size = rel(0.85), hjust = 1),
        plot.title = element_text(hjust=0.5)) +
  coord_fixed() +
  labs(x='Sample Type',y='Score',title='xCell')

#################################
# if plotting all at once:

# library(cowplot)
# plot_grid(g1,g2,g3,nrow=3)

# g1 <- ggplot(x,aes(reorder(x2$`Sample Type`,x2$value,mean,order=TRUE),variable,fill=value)) +
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
#                        name="Mean\nAbundance") +
#   theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 70, #vjust = 1,
#                                    size = rel(0.85), hjust = 1),
#         axis.text.y = element_text(size = rel(0.85)),
#         plot.title = element_text(hjust=0.5)) +
#   coord_fixed() +
#   labs(x='Sample Type',y='Score',title='CIBERSORT - Relative') +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# g3 <- ggplot(x3,aes(reorder(x2$`Sample Type`,x2$value,mean,order=TRUE),variable,fill=value)) +
#   geom_tile(color = "white")+
#   scale_fill_gradient(low = "white", high = "red",
#                        # midpoint = 0, limit = c(min(df.sub.melt$value),max(df.sub.melt$value)), space = "Lab",
#                        name="Mean\nAbundance") +
#   theme_minimal()+ # minimal theme
#   theme(axis.text.x = element_text(angle = 70, #vjust = 1,
#                                    size = rel(0.85), hjust = 1),
#         axis.text.y = element_text(size = rel(0.85)),
#         plot.title = element_text(hjust=0.5),
#         legend.title = element_text(size = rel(0.85))) +
#   coord_fixed() +
#   labs(x='Sample Type',y='Score',title='xCell')
# plot_grid(g1,g2,g3,nrow=3,rel_heights = 0.2)

# independently save g1, g2, g3
# g3
scale <- 3
tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','xcell','.png'),width=1000*scale,height=400*scale,res=100*scale,units="px")
print(g3)
dev.off()

tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','cibersort-rel','.png'),width=1000*scale,height=400*scale,res=100*scale,units="px")
print(g1)
dev.off()
tiff(paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/','cibersort-abs','.png'),width=1000*scale,height=400*scale,res=100*scale,units="px")
print(g2)
dev.off()


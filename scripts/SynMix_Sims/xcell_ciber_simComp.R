# Written by Andrew Marderstein (2018-2019). Contact: anm2868@med.cornell.edu

# Script for comparing xcell vs cibersort in simulation performance

library(data.table)

# on cluster
ciber.rel <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/SynMix_1.CIBERSORT.ABS-F.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
ciber.abs <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/SynMix_1.CIBERSORT.ABS-T.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
xcell <- fread('/athena/elementolab/scratch/anm2868/GTEx/gtexInfilSim/processed/XCell_simoutput_merged.txt',data.table = F,stringsAsFactors = F)

# local
ciber.rel <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/SynMix_Sims/SynMix_1.CIBERSORT.ABS-F.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
ciber.abs <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/SynMix_Sims/SynMix_1.CIBERSORT.ABS-T.QN-F.perm-1000.V2.txt',data.table = F,stringsAsFactors = F)
xcell <- fread('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/SynMix_Sims/XCell_simoutput_merged.txt',data.table = F,stringsAsFactors = F)

xcell.normalized=FALSE
if (xcell.normalized) {xcell[,2:65] <- xcell[,2:65]/apply(xcell[,2:65],1,sum)}
xcell$CD4_Tcells_Sum <- apply(xcell[,c('CD4+ memory T-cells',
                                       'CD4+ naive T-cells',
                                       'CD4+ T-cells',
                                       'CD4+ Tcm',
                                       'CD4+ Tem',
                                       'Th1 cells',
                                       'Th2 cells',
                                       'Tregs'
                                       )],1,sum)
xcell$CD8_Tcells_Sum <- apply(xcell[,c('CD8+ naive T-cells',
                                       'CD8+ T-cells',
                                       'CD8+ Tcm',
                                       'CD8+ Tem')],1,sum)
ciber.rel$CD4_Tcells_Sum <- apply(ciber.rel[,c('T cells CD4 naive',
                                           'T cells CD4 memory resting',
                                           'T cells CD4 memory activated',
                                           'T cells follicular helper',
                                           'T cells regulatory (Tregs)')],1,sum)
ciber.abs$CD4_Tcells_Sum <- apply(ciber.abs[,c('T cells CD4 naive',
                                           'T cells CD4 memory resting',
                                           'T cells CD4 memory activated',
                                           'T cells follicular helper',
                                           'T cells regulatory (Tregs)')],1,sum)
ciber.abs$CD8_Tcells_Sum <- ciber.abs[,'T cells CD8']
ciber.rel$CD8_Tcells_Sum <- ciber.rel[,'T cells CD8']

ciber.abs$`CD8_Tcells_Sum` <- ciber.abs$`T cells CD8`
ciber.rel$`CD8_Tcells_Sum` <- ciber.rel$`T cells CD8`

cor(xcell$`CD4_Tcells_Sum`,xcell$CD4_Naive_Actual_Pct)
cor(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual_Pct)
cor(xcell$`CD4_Tcells_Sum`,xcell$CD4_Naive_Actual)
cor(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual)

cor(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual_Pct)
cor(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual_Pct)
cor(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual)
cor(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual)

cor(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual_Pct)
cor(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual_Pct)
cor(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual)
cor(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual)

cor(ciber.abs$`CD8_Tcells_Sum`,xcell$`CD8_Tcells_Sum`)
cor(ciber.abs$`CD4_Tcells_Sum`,xcell$`CD4_Tcells_Sum`)


# for plotting/visualization purposes
# plot(xcell$`CD4_Tcells_Sum`,xcell$CD4_Naive_Actual_Pct)
# plot(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual_Pct)
# plot(xcell$CD4_Naive_Actual,xcell$`CD4_Tcells_Sum`)
# plot(xcell$`CD8_Tcells_Sum`,xcell$CD8_Naive_Actual)
# plot(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual_Pct)
# plot(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual_Pct)
# plot(ciber.abs$`CD4_Tcells_Sum`,ciber.abs$CD4_Naive_Actual)
# plot(ciber.abs$`CD8_Tcells_Sum`,ciber.abs$CD8_Naive_Actual)
# plot(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual_Pct)
# plot(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual_Pct)
# plot(ciber.rel$`CD4_Tcells_Sum`,ciber.rel$CD4_Naive_Actual)
# plot(ciber.rel$`CD8_Tcells_Sum`,ciber.rel$CD8_Naive_Actual)
# plot(ciber.abs$`CD8_Tcells_Sum`,xcell$`CD8_Tcells_Sum`)
# plot(ciber.abs$`CD4_Tcells_Sum`,xcell$`CD4_Tcells_Sum`)

# case example, cibersort vs xcell:
df <- data.frame(actual=xcell$CD4_Naive_Actual,xcell=xcell$`CD4_Tcells_Sum`,ciber.rel=ciber.rel$`CD4_Tcells_Sum`,ciber.abs=ciber.abs$`CD4_Tcells_Sum`)
ggplot(df,aes(x=actual,y=xcell,col=ciber.abs)) + geom_point()
x <- subset(df,actual==10)
x$ID <- 1:nrow(x)
x <- x[c(1,7),]
x[2,] <- x[2,]/x[1,]
x[1,] <- x[1,]/x[1,]
x <- subset(reshape2::melt(x,id.vars='ID'),variable %in% c('xcell','ciber.abs'))
ggplot(x,aes(x=as.factor(variable),group=as.factor(ID),fill=as.factor(ID),y=value)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title='synthetic mixes w/ equal immune content',x='Method',y='Relative score') +
  scale_fill_manual(values=c('steelblue2','orange2')) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5)) + theme(legend.position='none') +
  scale_x_discrete(labels=c('xCell','CIBERSORT-Absolute'))
xcell[c(1,58),]

# snp analysis simulations:
df <- data.frame(xcell=xcell$`CD4_Tcells_Sum`,ciber.rel=ciber.rel$`CD4_Tcells_Sum`,ciber.abs=ciber.abs$`CD4_Tcells_Sum`)
nsim <- 10000
p.df <- matrix(NA,nsim,4)
Y.df <- matrix(NA,nrow(df),3)
no_eff <- TRUE #FALSE
for (j in 1:nsim) {
  G <- rbinom(nrow(df),2,0.4)
  p <- rep(NA,3)
  var.g.vec <- runif(3,0,0.1)
  for (i in 1:3) {
    if (!no_eff) {
      var.g <- var.g.vec[i]
      Y.G <- G
      Y.env <- df[,i]
      Y.env <- Y.env * sqrt(var(Y.G) / var(Y.env))
      Y.G <- Y.G * sqrt(var.g)
      Y <- Y.G + Y.env
    } else if (no_eff) {
      Y.env <- df[,i]
      Y <- Y.env
    }
    Y.df[,i] <- Y
    mod <- (lm(Y ~ G))
    p[i] <- summary(mod)$coef[2,4]
    
  }
  df.pheno2 <- as.data.frame(t(Y.df))
  cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno2)
  print('Running sped up pre-calculated covariance matrix method...')
  p.df[j,] <- c(p,
        empiricalBrownsMethod2(data_matrix = Y.df,p_values = p,extra_info = F,covar_matrix = cov.matrix)
  )
}
p.df <- as.data.frame(p.df)
colnames(p.df) <- c('xcell','cib.rel','cib.abs','ebm')
for (i in 1:4) {print(sum(p.df[,i] < 0.05))}
for (i in 1:4) {print(sum(p.df[,i] < 0.05/10000))}
for (i in 1:4) {print(sum(p.df[,i] < 5e-8))}

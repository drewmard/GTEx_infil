# initialize:
id <- list(); 
id[[1]] <- 166
id[[2]] <- 56
id[[3]] <- 81

# 3: plot genotype phenotype plots for top 3 hits
f <- '/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/gtex_all.filter.name.txt'
df.pheno <- fread(f,stringsAsFactors = F,data.table = F)

snp <- list(); rs <- list(); g.geno <- list(); phenoName <- list(); df.snp <- list(); pheno <- list(); df.geno <- list()
tis <- list('Thyroid','Colon - Sigmoid','Esophagus - Muscularis')
phenoName <- list("Tfh cells (trans resid)",
                  "Lymphocytes (trans resid)",
                  "Mast cells (trans resid)")
pheno <- as.list(paste0('pheno',id))

snp <- list('10_22337752_A_G_b37','22_23961126_C_T_b37','17_77995143_A_G_b37')
rs <- list('rs6482199','rs56234965','rs9989443')
A1 <- list(); A2 <- list()

i=1
f <- paste0('/Users/andrewmarderstein/Documents/Research/GTEx/Infiltration/GTEx_infil/output/GeneticAnalysis/Genotype_Subset/gtex_all.filter.name.',snp[[i]],'.raw')
df.snp[[i]] <- fread(f,data.table = F,stringsAsFactors = F)
df.geno[[i]] <- merge(df.pheno[,c('IID',paste0(pheno[[i]],'.',1:3))],df.snp[[i]][,c(2,grep(snp[[i]],colnames(df.snp[[i]])))],by='IID')
df.geno[[i]] <- subset(df.geno[[i]],!(is.na(df.geno[[i]][,grep(snp[[i]],colnames(df.geno[[i]]))]) | (df.geno[[i]][,paste0(pheno[[i]],'.',1)] == -9)))

plot(df.geno[[i]][,'pheno166.1'],df.geno[[i]][,'pheno166.2'])
plot(df.geno[[i]][,'pheno166.2'],df.geno[[i]][,'pheno166.3'])
ggplot(df.geno[[i]],aes(x=pheno166.2,y=pheno166.3)) + geom_point() + geom_smooth()

cor(df.geno[[i]][,-1])



# actual
i=1
df <- df.geno[[i]]
Y.df <- df[,2:(ncol(df)-1)]
p <- rep(NA,3)
for (k in 1:3) {
  mod <- (lm(Y.df[,k] ~ G))
  p[k] <- summary(mod)$coef[2,4]
}
df.pheno2 <- as.data.frame(t(Y.df))
cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno2)
# print('Running sped up pre-calculated covariance matrix method...')
p.real <- c(p,
              empiricalBrownsMethod2(data_matrix = Y.df,p_values = p,extra_info = F,covar_matrix = cov.matrix)
)

nsim=10000
p.df <- matrix(NA,nsim,4)
G <- df[,ncol(df)]
for (j in 1:nsim) {
  
  if (j %%100 == 0) {print(j)}
  i <- sample(1:nrow(df),size = nrow(df),replace = F)
  df[,2:(ncol(df)-1)] <- df[i,2:(ncol(df)-1)] 
  
  
  Y.df <- df[,2:(ncol(df)-1)]
  # cov(Y.df)
  p <- rep(NA,3)
  for (k in 1:3) {
    mod <- (lm(Y.df[,k] ~ G))
    p[k] <- summary(mod)$coef[2,4]
  }
  df.pheno2 <- as.data.frame(t(Y.df))
  cov.matrix <- EmpiricalBrownsMethod:::calculateCovariances(df.pheno2)
  # print('Running sped up pre-calculated covariance matrix method...')
  p.df[j,] <- c(p,
                empiricalBrownsMethod2(data_matrix = Y.df,p_values = p,extra_info = F,covar_matrix = cov.matrix)
  )
}
p.df <- as.data.frame(p.df)
colnames(p.df) <- c('cib.rel','cib.abs','xcell','ebm')
for (i in 1:4) {print(sum(p.df[,i] < 0.05))}
for (i in 1:4) {print(sum(p.df[,i] < 0.05/10000))}
for (i in 1:4) {print(sum(p.df[,i] < 5e-8))}

p.df[c(which.min(p.df[,1]),which.min(p.df[,2]),which.min(p.df[,3]),which.min(p.df[,4])),]




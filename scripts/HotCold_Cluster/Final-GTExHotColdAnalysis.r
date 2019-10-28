# Written by Manik Uppal (2018-2019). Contact: mdu2002@med.cornell.edu

library("openxlsx")
library("pheatmap")
library("RColorBrewer")
library("gplots")
library("edgeR")
library("limma")
library("DESeq2")
library("ggplot2")

###################################################################################
###################################################################################
###################################################################################
#BEGIN NORMAL SCRIPT

path = "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/eQTL/RawData/"

CountsRaw = paste(path,"All_Tissue_Site_Details.combined.reads.txt",sep="")
CountsRaw_mod = paste(path,"GTEx-Counts_AvgTfd.txt",sep="")
TPMsRaw_mod = paste(path,"GTEx-TPMs_AvgTfd.txt",sep="")
xCellCounts = paste(path,"xCellOutput-GTEx-TPM.txt",sep="")
annoCounts = paste(path,"GTExSampleAnno_mod.csv",sep="")
totalAnno = paste(path,"GTExAnno_all.xlsx",sep="")

raw_dat <- read.table(CountsRaw_mod,header=TRUE)
row.names <- raw_dat$Gene
rownames(raw_dat) <- row.names
raw_dat <- raw_dat[,-1]

cellTypes.df <- data.frame(ciber=c('T cells CD8','T cells CD4 naive','CD4_memory','Neutrophils','MacrophageSum','Bcellsum','NK_Sum','DendriticSum','MastSum','Myeloid_Sum','T cells follicular helper','T cells regulatory (Tregs)','T cells gamma delta','Monocytes','Eosinophils','Lymph_Sum'), xcell=c('CD8Sum','CD4+ naive T-cells','CD4_memory','Neutrophils','MacrophageSum','Bcellsum','NK cells','DendriticSum','Mast cells','Myeloid_Sum','Th_Sum','Tregs','Tgd cells','Monocytes','Eosinophils','Lymph_Sum'),stringsAsFactors = F)

print("Loading xCell Deconvolution data...")
xCell.orig <- read.csv("/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/MiscWorkup/XCell.all_tissues.txt",stringsAsFactors=FALSE,sep="\t",header=F)
xCell.orig[,69] <- gsub('-','.',xCell.orig[,69])
xCell <- xCell.orig
colnames(xCell) <- xCell[1,]
xCell <- xCell[-1,c(2:69,74:83)]
xCell$SAMP <- gsub('-','.',xCell$SAMP)
xCellSamples <- xCell$SAMP
xCell <- xCell[,-68]
xCell <- sapply(xCell, as.numeric)
rownames(xCell) <- xCellSamples
deconv_data.xCell <- t(xCell)

print("Loading Cibersort Deconvolution data...")
cibersortAbs.orig <- read.csv("/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/MiscWorkup/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-T.QN-F.perm-1000.txt",stringsAsFactors=FALSE,sep="\t")
cibersortAbs <- cibersortAbs.orig
cibersortAbs <- cibersortAbs[,-c(1,25:32)]
cibersortAbs$Input.Sample <- gsub('-','.',cibersortAbs$Input.Sample)
rownames(cibersortAbs) <- cibersortAbs$Input.Sample
cibersortAbs <- cibersortAbs[,-1]
deconv_data.Cibersort <- t(cibersortAbs)

cibersortRel.orig <- read.csv("/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/MiscWorkup/GTEx_v7_genexpr_ALL.CIBERSORT.ABS-F.QN-F.perm-1000.txt",stringsAsFactors=FALSE,sep="\t")
cibersortRel <- cibersortRel.orig
cibersortRel <- cibersortRel[,-c(1,25:32)]
cibersortRel$Input.Sample <- gsub('-','.',cibersortRel$Input.Sample)
rownames(cibersortRel) <- cibersortRel$Input.Sample
cibersortRel <- cibersortRel[,-1]
deconv_data.CibersortRel <- t(cibersortRel)

print("Loading Annotation data...")
total_anno <- read.xlsx(totalAnno)
total_anno$SAMPID <- gsub('-','.',total_anno$SAMPID)
rownames(total_anno) <- total_anno$SAMPID
total_anno <- total_anno[,-1]

print("Matching sample order of Annotations and Deconvolution data...")
annoInd <- which(rownames(total_anno) %in% colnames(deconv_data.xCell))
total_anno_mod <- total_anno[annoInd,]
seIdx <- match(colnames(deconv_data.xCell), rownames(total_anno_mod))
total_anno_mod <- total_anno_mod[seIdx,]
seIdx <- match(colnames(deconv_data.xCell), colnames(deconv_data.Cibersort))
deconv_data.Cibersort <- deconv_data.Cibersort[,seIdx]
seIdx <- match(colnames(deconv_data.xCell), colnames(deconv_data.CibersortRel))
deconv_data.CibersortRel <- deconv_data.CibersortRel[,seIdx]


colsToRemove <- c("Adipocytes","Astrocytes","Chondrocytes","CLP","CMP","Endothelial cells","Epithelial cells","Erythrocytes","Fibroblasts","GMP","Hepatocytes","HSC","Keratinocytes","ly Endothelial cells","Megakaryocytes","Melanocytes","MEP","Mesangial cells","MPP","MSC","mv Endothelial cells","Myocytes","Neurons","Osteoblast","Pericytes","Platelets","Preadipocytes","Sebocytes","Skeletal muscle","Smooth muscle","StromaScore","ImmuneScore","MicroenvironmentScore")
indices <- which(rownames(deconv_data.xCell) %in% colsToRemove)
upd_deconv.xCell <- deconv_data.xCell[-indices, ]
deconv_data.xCell <- upd_deconv.xCell
#This function returns a dataframe with the values (FPKM, enrichment, etc.) on organ-level basis; organ must be specified
isolateOrganComparison <- function(sampleVals, anno, OrganString)
{
    #since I generated a matching ordered deconv/raw data and annotation data with matchSampleOrder, let's just 1) isolate the subset of samples with Tissue = OrganString in "anno", 2) generate the vector of names of samples, 3) isolate the subset of deconv/raw data pertaining to those sampleIDs, 4) return those sampleIDs
    anno_subset <- anno[which(anno$SMTSD == OrganString),]
    colIndices <- which(colnames(sampleVals) %in% rownames(anno_subset))
    sampleVals_subset <- sampleVals[,colIndices]
    return(sampleVals_subset)
    
}

Tissues <- unique(total_anno_mod$SMTSD)
Tissues <- Tissues[-54]
anno_dat <- total_anno_mod

###################################################################################
###################################################################################
###################################################################################
#Consensus Clustering workflow
#The following loads the 189 tissue-cell pairs that are analyzed
print("Loading Tissue Cell Pairs")
pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"
tissCellPath <- paste(pureCellFilePath,"TissueCellPairs-Updated.txt",sep="")
tissCellPairs <- read.table(tissCellPath, header=T, sep="\t",stringsAsFactors=FALSE)

#Consensus Clustering using triplicates
for(p in 1:length(tissCellPairs$Tissue))
{
    TissCell <- tissCellPairs[p,]
    OrganString <- as.character(TissCell[1])
    ImmCellType.Cibersort <- as.character(TissCell[2])
    ImmCellType.xCell <- as.character(TissCell[3])
    print(OrganString)
    print(ImmCellType.xCell)
    
    #xCell
    organ_deconv.xCell <- isolateOrganComparison(deconv_data.xCell,anno_dat,OrganString)
    rowIndImm.xCell <- which(rownames(organ_deconv.xCell) %in% ImmCellType.xCell)
    ImmCell.organ.xCell <- as.vector(organ_deconv.xCell[rowIndImm.xCell,])
    attr(ImmCell.organ.xCell,"names") <- colnames(organ_deconv.xCell)
    triplicate.xCell <- rbind(ImmCell.organ.xCell,ImmCell.organ.xCell,ImmCell.organ.xCell)
    ConsClustFileName <- paste("PureCell-",OrganString,"-",ImmCellType.xCell,"-xCell",sep="")
    #print("xCell")
    result.xCell <- ConsensusCluster(as.matrix(triplicate.xCell), maxK = 20, reps = 2000,clusterAlg = "km",title = ConsClustFileName,verbose=FALSE,plot="pngBMP",writeTable=TRUE)
    #print(dim(triplicate.xCell))
    
    #Cibersort Absolute
    organ_deconv.Cibersort <- isolateOrganComparison(deconv_data.Cibersort,anno_dat,OrganString)
    rows.Cibersort <- rownames(organ_deconv.Cibersort)
    rows.Cibersort <- gsub("\\.", " ",rows.Cibersort)
    rowIndImm.Cibersort <- which(rows.Cibersort %in% ImmCellType.Cibersort)
    ImmCell.organ.Cibersort <- as.vector(organ_deconv.Cibersort[rowIndImm.Cibersort,])
    attr(ImmCell.organ.Cibersort,"names") <- colnames(organ_deconv.Cibersort)
    triplicate.Cibersort <- rbind(ImmCell.organ.Cibersort,ImmCell.organ.Cibersort,ImmCell.organ.Cibersort)
    ConsClustFileName <- paste("PureCell-",OrganString,"-",ImmCellType.Cibersort,"-Cibersort",sep="")
    #print("Cibersort")
    result.Cibersort <- ConsensusCluster(as.matrix(triplicate.Cibersort), maxK = 20, reps = 2000,clusterAlg = "km",title = ConsClustFileName,verbose=FALSE,plot="pngBMP",writeTable=TRUE)
    #print(dim(triplicate.Cibersort))
    
    #Cibersort Relative
    organ_deconv.Cibersort <- isolateOrganComparison(deconv_data.CibersortRel,anno_dat,OrganString)
    rows.Cibersort <- rownames(organ_deconv.Cibersort)
    rows.Cibersort <- gsub("\\.", " ",rows.Cibersort)
    rowIndImm.Cibersort <- which(rows.Cibersort %in% ImmCellType.Cibersort)
    ImmCell.organ.Cibersort <- as.vector(organ_deconv.Cibersort[rowIndImm.Cibersort,])
    attr(ImmCell.organ.Cibersort,"names") <- colnames(organ_deconv.Cibersort)
    triplicate.Cibersort <- rbind(ImmCell.organ.Cibersort,ImmCell.organ.Cibersort,ImmCell.organ.Cibersort)
    ConsClustFileName <- paste("PureCell-",OrganString,"-",ImmCellType.Cibersort,"-CibersortRel",sep="")
    #print("Cibersort")
    result.Cibersort <- ConsensusCluster(as.matrix(triplicate.Cibersort), maxK = 20, reps = 2000,clusterAlg = "km",title = ConsClustFileName,verbose=FALSE,plot="pngBMP",writeTable=TRUE)
}

###################################################################################
###################################################################################
###################################################################################
pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"

tissCellPath <- paste(pureCellFilePath,"TissueCellPairs-Updated.txt",sep="")
tissCellPairs <- read.table(tissCellPath, header=T, sep="\t",stringsAsFactors=FALSE)
k_mat <- read.table(paste(pureCellFilePath,"TissCellPairs-k_opt-Updated.txt",sep=""),header=TRUE,stringsAsFactors=FALSE,sep="\t")
k_opt.xCell <- k_mat$`k.xCell`
k_opt.Cibersort <- k_mat$`k.Cibersort.Abs`
k_opt.CibersortRel <- k_mat$`k.Cibersort.Rel`
hotcoldMat <- matrix(, nrow = length(tissCellPairs$Tissue), ncol = 2)

#Produce Hot Cold Matrices with Full Covariate Assignments
for(p in 1:length(tissCellPairs$Tissue))
{
    TissCell <- tissCellPairs[p,]
    OrganString <- as.character(TissCell[1])
    ImmCellType.Cibersort <- as.character(TissCell[2])
    ImmCellType.xCell <- as.character(TissCell[3])
    ImmCellType.CibersortRel <-ImmCellType.Cibersort
    
    #load xCell hot and cold results
    organ_deconv.xCell <- isolateOrganComparison(deconv_data.xCell,anno_dat,OrganString)
    rowIndImm.xCell <- which(rownames(organ_deconv.xCell) %in% ImmCellType.xCell)
    ImmCell.organ.xCell <- as.vector(organ_deconv.xCell[rowIndImm.xCell,])
    attr(ImmCell.organ.xCell,"names") <- colnames(organ_deconv.xCell)
    k.xCell <- k_opt.xCell[p]
    
    ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.xCell,"-xCell",sep="")
    consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.xCell,".consensusClass.csv",sep="")
    optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
    colnames(optClusters) <- c("Sample","Cluster")
    rownames(optClusters) <- optClusters$Sample
    
    tree <- optClusters
    numClusters = length(unique(tree$Cluster))
    geneClustAvgs <- numeric(numClusters)
    for(l in 1:numClusters)
    {
        clustSamples <- as.vector(tree[tree$Cluster == l,]$Sample)
        clust_SentGeneExpr <- ImmCell.organ.xCell[which(labels(ImmCell.organ.xCell) %in% clustSamples)]
        clust_SentGeneExpr.avg <- mean(as.numeric(clust_SentGeneExpr))
        geneClustAvgs[l] <- clust_SentGeneExpr.avg
    }
    nonInflBranch <- which.min(geneClustAvgs)
    InflBranch <- which.max(geneClustAvgs)
    branches.keep <- c(nonInflBranch, InflBranch)
    tree.subset.xCell <- tree[which(tree$Cluster %in% branches.keep),]
    numInfl <- nrow(tree.subset.xCell[tree.subset.xCell$Cluster == InflBranch,])
    numNonInfl <- nrow(tree.subset.xCell[tree.subset.xCell$Cluster == nonInflBranch,])
    
    for(q in 1:nrow(tree.subset.xCell))
    {
        sample <- as.character(tree.subset.xCell[q,1])
        sampVec <- unlist(strsplit(as.character(sample), "[.]"))
        newName <- paste(sampVec[1],sampVec[2],sep="-")
        tree[q,1] <- newName
        cluster <- tree.subset.xCell[q,2]
        if(cluster == InflBranch)
        {
            tree.subset.xCell[q,2] <- 1
        }
        else if(cluster != InflBranch)
        {
            tree.subset.xCell[q,2] <- 0
        }
    }
    
    ###########################################
    ###########################################
    ###########################################
    
    
    ####CibersortAbsolute
    organ_deconv.Cibersort <- isolateOrganComparison(deconv_data.Cibersort,anno_dat,OrganString)
    k.Cibersort <- k_opt.Cibersort[p]
    rows.Cibersort <- rownames(organ_deconv.Cibersort)
    rows.Cibersort <- gsub("\\.", " ",rows.Cibersort)
    rowIndImm.Cibersort <- which(rows.Cibersort  %in% ImmCellType.Cibersort)
    ImmCell.organ.Cibersort <- as.vector(organ_deconv.Cibersort[rowIndImm.Cibersort,])
    attr(ImmCell.organ.Cibersort,"names") <- names(organ_deconv.Cibersort[rowIndImm.Cibersort,])
    
    ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.Cibersort,"-CibersortAbs",sep="")
    consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.Cibersort,".consensusClass.csv",sep="")
    optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
    colnames(optClusters) <- c("Sample","Cluster")
    rownames(optClusters) <- optClusters$Sample
    
    tree <- optClusters
    numClusters = length(unique(tree$Cluster))
    geneClustAvgs <- numeric(numClusters)
    for(l in 1:numClusters)
    {
        clustSamples <- as.vector(tree[tree$Cluster == l,]$Sample)
        clust_SentGeneExpr <- ImmCell.organ.Cibersort[which(labels(ImmCell.organ.Cibersort) %in% clustSamples)]
        clust_SentGeneExpr.avg <- mean(as.numeric(clust_SentGeneExpr))
        geneClustAvgs[l] <- clust_SentGeneExpr.avg
    }
    nonInflBranch <- which.min(geneClustAvgs)
    InflBranch <- which.max(geneClustAvgs)
    branches.keep <- c(nonInflBranch, InflBranch)
    tree.subset.Cibersort <- tree[which(tree$Cluster %in% branches.keep),]
    numInfl <- nrow(tree.subset.Cibersort[tree.subset.Cibersort$Cluster == InflBranch,])
    numNonInfl <- nrow(tree.subset.Cibersort[tree.subset.Cibersort$Cluster == nonInflBranch,])
    
    for(q in 1:nrow(tree.subset.Cibersort))
    {
        sample <- as.character(tree.subset.Cibersort[q,1])
        sampVec <- unlist(strsplit(as.character(sample), "[.]"))
        newName <- paste(sampVec[1],sampVec[2],sep="-")
        tree[q,1] <- newName
        cluster <- tree.subset.Cibersort[q,2]
        if(cluster == InflBranch)
        {
            tree.subset.Cibersort[q,2] <- 1
        }
        else if(cluster != InflBranch)
        {
            tree.subset.Cibersort[q,2] <- 0
        }
    }
    
    ###########################################
    ###########################################
    ###########################################
    
    ####CibersortRelative
    organ_deconv.CibersortRel <- isolateOrganComparison(deconv_data.CibersortRel,anno_dat,OrganString)
    k.CibersortRel <- k_opt.CibersortRel[p]
    rows.CibersortRel <- rownames(organ_deconv.CibersortRel)
    rows.CibersortRel <- gsub("\\.", " ",rows.CibersortRel)
    rowIndImm.CibersortRel <- which(rows.CibersortRel  %in% ImmCellType.CibersortRel)
    ImmCell.organ.CibersortRel <- as.vector(organ_deconv.CibersortRel[rowIndImm.CibersortRel,])
    attr(ImmCell.organ.CibersortRel,"names") <- names(organ_deconv.CibersortRel[rowIndImm.CibersortRel,])
    
    ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.CibersortRel,"-CibersortRel",sep="")
    consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.CibersortRel,".consensusClass.csv",sep="")
    optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
    colnames(optClusters) <- c("Sample","Cluster")
    rownames(optClusters) <- optClusters$Sample
    
    tree <- optClusters
    numClusters = length(unique(tree$Cluster))
    geneClustAvgs <- numeric(numClusters)
    for(l in 1:numClusters)
    {
        clustSamples <- as.vector(tree[tree$Cluster == l,]$Sample)
        clust_SentGeneExpr <- ImmCell.organ.CibersortRel[which(labels(ImmCell.organ.CibersortRel) %in% clustSamples)]
        clust_SentGeneExpr.avg <- mean(as.numeric(clust_SentGeneExpr))
        geneClustAvgs[l] <- clust_SentGeneExpr.avg
    }
    nonInflBranch <- which.min(geneClustAvgs)
    InflBranch <- which.max(geneClustAvgs)
    branches.keep <- c(nonInflBranch, InflBranch)
    tree.subset.CibersortRel <- tree[which(tree$Cluster %in% branches.keep),]
    numInfl <- nrow(tree.subset.CibersortRel[tree.subset.CibersortRel$Cluster == InflBranch,])
    numNonInfl <- nrow(tree.subset.CibersortRel[tree.subset.CibersortRel$Cluster == nonInflBranch,])
    
    for(q in 1:nrow(tree.subset.CibersortRel))
    {
        sample <- as.character(tree.subset.CibersortRel[q,1])
        sampVec <- unlist(strsplit(as.character(sample), "[.]"))
        newName <- paste(sampVec[1],sampVec[2],sep="-")
        tree[q,1] <- newName
        cluster <- tree.subset.CibersortRel[q,2]
        if(cluster == InflBranch)
        {
            tree.subset.CibersortRel[q,2] <- 1
        }
        else if(cluster != InflBranch)
        {
            tree.subset.CibersortRel[q,2] <- 0
        }
    }
    
    ###########################################
    ###########################################
    ###########################################
    
    hot.xCell <- tree.subset.xCell[which(tree.subset.xCell$Cluster == 1),1]
    hot.Cibersort <- tree.subset.Cibersort[which(tree.subset.Cibersort$Cluster == 1),1]
    hot.CibersortRel <- tree.subset.CibersortRel[which(tree.subset.CibersortRel$Cluster == 1),1]
    consHotSamps <- intersect(intersect(hot.xCell,hot.Cibersort),hot.CibersortRel)
    cold.xCell <- tree.subset.xCell[which(tree.subset.xCell$Cluster == 0),1]
    cold.Cibersort <- tree.subset.Cibersort[which(tree.subset.Cibersort$Cluster == 0),1]
    cold.CibersortRel <- tree.subset.CibersortRel[which(tree.subset.CibersortRel$Cluster == 0),1]
    consColdSamps <- intersect(intersect(cold.xCell,cold.Cibersort),cold.CibersortRel)
    
    consSamples <- c(consHotSamps,consColdSamps)
    tree.final <-tree.subset.xCell[which(tree.subset.xCell$Sample %in% consSamples),]
    heatString <- character(0)
    for(q in 1:nrow(tree.final))
    {
        if(tree.final[q,2] == 0)
        {
            heatString <- c(heatString,"Cold")
        }
        if(tree.final[q,2] == 1)
        {
            heatString <- c(heatString,"Hot")
        }
    }
    tree.final$Infiltration <- heatString
    anno_inds <- which(rownames(anno_dat) %in% tree.final$Sample)
    autolysisScore <- anno_dat[anno_inds,]$SMATSSCR
    attr(autolysisScore,"names") <- rownames(anno_dat)[anno_inds]
    seIdx <- match(tree.final$Sample, names(autolysisScore))
    autolysisScore <- autolysisScore[seIdx]
    tree.final$Autolysis <- autolysisScore
    
    collectionSite <- anno_dat[anno_inds,]$SMCENTER
    attr(collectionSite,"names") <- rownames(anno_dat)[anno_inds]
    seIdx <- match(tree.final$Sample, names(collectionSite))
    collectionSite <- collectionSite[seIdx]
    tree.final$Site <- collectionSite
    
    anno_inds <- which(xCell.orig[,69] %in% tree.final$Sample)
    covariates <- xCell.orig[anno_inds,c(69,71:73)]
    seIdx <- match(tree.final$Sample, covariates[,1])
    covariates <- covariates[seIdx,]
    tree.final$SEX <- covariates[,2]
    tree.final$AGE <- covariates[,3]
    tree.final$DTHHRDY <- covariates[,4]
    tree.subset <- tree.final
    write.table(tree.subset,file = paste(pureCellFilePath,"HotColdMatrixModels/PureCell-",OrganString,"-",ImmCellType.xCell,"-Consensus","-HotColdMatrix.txt",sep=""),quote = FALSE,row.names=FALSE,sep="\t")
    
    hotcoldMat[p,1] <- table(tree.final$Cluster)[1]
    hotcoldMat[p,2] <- table(tree.final$Cluster)[2]
    
    
    
}
ConsensusHotColdNumbers <- data.frame(Tissue=tissCellPairs$Tissue,CellType=tissCellPairs$xCell,Cold=hotcoldMat[,1],Hot=hotcoldMat[,2])
write.table(ConsensusHotColdNumbers,paste(pureCellFilePath,"ConsensusHotColdNumbers-New.txt",sep=""),sep="\t",row.names=F,quote=F)

#Create Covariate Models and Perform DGE
tissCellModel <- character(0)
sampleClustersAndCovariates <- character(0)
DeconvDGEPath <- paste(pureCellFilePath,"PureCellConsensusDGE/",sep="")
#DeconvDGEPath <- "/Users/ManikUppal/Desktop/ValueRegressed/"
for(p in 1:length(tissCellPairs$Tissue))
{
    TissCell <- tissCellPairs[p,]
    OrganString <- as.character(TissCell[1])
    ImmCellType.Cibersort <- as.character(TissCell[2])
    ImmCellType.xCell <- as.character(TissCell[3])
    ImmCellType.CibersortRel <-ImmCellType.Cibersort
    organ_dat <- isolateOrganComparison(raw_dat,anno_dat,OrganString)
    organ_dat.rd <- as.matrix(as.data.frame(lapply(organ_dat, as.integer)))
    rownames(organ_dat.rd) <- rownames(organ_dat)
    NumSamples <- ncol(organ_dat)
    
    tree.subset <- read.table(paste(pureCellFilePath,"HotColdMatrixModels/PureCell-",OrganString,"-",ImmCellType.xCell,"-Consensus","-HotColdMatrix.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")
    samples.subset <- as.vector(tree.subset$Sample)
    geneClustComparison.df <- organ_dat.rd[,which(colnames(organ_dat.rd) %in% samples.subset)]
    geneClustComparison.df <- geneClustComparison.df[apply(geneClustComparison.df, MARGIN = 1, FUN = function(x) sd(x) != 0),]
    seIdx <- match(colnames(geneClustComparison.df), tree.subset$Sample)
    tree.subset <- tree.subset[seIdx, ]
    rownames(tree.subset) <- tree.subset$Sample
    tree.subset$Infiltration <- factor(tree.subset$Infiltration)
    tree.subset$Infiltration <- relevel(tree.subset$Infiltration,"Cold")
    tree.subset$AGE <- as.numeric(as.factor(tree.subset$AGE))#
    tree.subset$DTHHRDY <- as.factor(tree.subset$DTHHRDY)
    tree.subset$SEX <- as.factor(tree.subset$SEX)
    tree.subset$Autolysis <- as.numeric(as.character(tree.subset$Autolysis))#
    tree.subset$Site <- as.factor(tree.subset$Site)
    
    if(p >= 185)
    {
        tree.subset <- tree.subset[,-4]
    }
    tree.subset <- tree.subset[rowSums(is.na(tree.subset)) == 0,]
    geneClustComparison.df <- geneClustComparison.df[,which(colnames(geneClustComparison.df) %in% tree.subset$Sample)]
    numNonInfl <- table(tree.subset$Cluster)[1]
    numInfl <- table(tree.subset$Cluster)[2]
    
    if(!is.na(numInfl)){
        if(numInfl >= 6 && numNonInfl >= 6)
        {
            MINTHRESHOLD <- 2
            HCDIFFTHRESHOLD <- 5
            modelString <- "tree.subset$Infiltration"
            
            
            sex <- table(tree.subset$Infiltration, tree.subset$SEX)
            if(p<185){
                autolysis <- table(tree.subset$Infiltration, tree.subset$Autolysis)}
            COD <- table(tree.subset$Infiltration, tree.subset$DTHHRDY)
            age <- table(tree.subset$Infiltration, tree.subset$AGE)
            site <- table(tree.subset$Infiltration, tree.subset$Site)
            
            ############
            #Need a minimum number of samples in each level to be included as covariate
            #SEX
            if(ncol(sex) ==2){
                if(length(which((sex[2,] > MINTHRESHOLD) == FALSE)) == 0)
                {
                    modelString <- paste("tree.subset$SEX + ",modelString,sep="")
                }
                if(length(which((sex[2,] > MINTHRESHOLD) == FALSE)) == 1)
                {
                    ind <- which(sex[2,] <= MINTHRESHOLD)
                    if(length(ind) != 0){
                        if( (sex[1,ind] - sex[2,ind] ) <= HCDIFFTHRESHOLD) {
                            modelString <- paste("tree.subset$SEX + ",modelString,sep="")
                        }}
                }}
            
            #Autolysis
            if(p<185){
                if(length(which((autolysis[2,] > MINTHRESHOLD) == FALSE)) == 0)
                {
                    modelString <- paste("tree.subset$Autolysis + ",modelString,sep="")
                }
                if(length(which((autolysis[2,] > MINTHRESHOLD) == FALSE)) == 1)
                {
                    ind <- which(autolysis[2,] <= MINTHRESHOLD)
                    if(length(ind) != 0){
                        if( (autolysis[1,ind] - autolysis[2,ind] ) <= HCDIFFTHRESHOLD) {
                            modelString <- paste("tree.subset$Autolysis + ",modelString,sep="")
                        }}
                }}
            
            #COD
            if(length(which((COD[2,] > MINTHRESHOLD) == FALSE)) == 0)
            {
                modelString <- paste("tree.subset$DTHHRDY + ",modelString,sep="")
            }
            if(length(which((COD[2,] > MINTHRESHOLD) == FALSE)) == 1)
            {
                ind <- which(COD[2,] <= MINTHRESHOLD)
                if(length(ind) != 0){
                    if( (COD[1,ind] - COD[2,ind] ) <= HCDIFFTHRESHOLD) {
                        modelString <- paste("tree.subset$DTHHRDY + ",modelString,sep="")
                    }}
            }
            
            #AGE
            if(length(which((age[2,] > MINTHRESHOLD) == FALSE)) == 0)
            {
                modelString <- paste("tree.subset$AGE + ",modelString,sep="")
            }
            if(length(which((age[2,] > MINTHRESHOLD) == FALSE)) == 1)
            {
                ind <- which(age[2,] <= MINTHRESHOLD)
                if(length(ind) != 0){
                    if( (age[1,ind] - age[2,ind] ) <= HCDIFFTHRESHOLD) {
                        modelString <- paste("tree.subset$AGE + ",modelString,sep="")
                    }}
            }
            
            #SITE
            if(length(which((site[2,] > MINTHRESHOLD) == FALSE)) == 0)
            {
                modelString <- paste("tree.subset$Site + ",modelString,sep="")
            }
            if(length(which((site[2,] > MINTHRESHOLD) == FALSE)) == 1)
            {
                ind <- which(site[2,] <= MINTHRESHOLD)
                if(length(ind) != 0){
                    if( (site[1,ind] - site[2,ind] ) <= HCDIFFTHRESHOLD) {
                        modelString <- paste("tree.subset$Site + ",modelString,sep="")
                    }}
            }
            
            ###########
            #comment out this block, it's only for testing using the xCell infiltration score as a covariate to see if tissue specific drivers can be uncovered
            #organ_deconv.xCell <- isolateOrganComparison(deconv_data.xCell,anno_dat,OrganString)
            #rowIndImm.xCell <- which(rownames(organ_deconv.xCell) %in% ImmCellType.xCell)
            #ImmCell.organ.xCell <- as.vector(organ_deconv.xCell[rowIndImm.xCell,])
            #attr(ImmCell.organ.xCell,"names") <- colnames(organ_deconv.xCell)
            
            #values<- ImmCell.organ.xCell[which(names(ImmCell.organ.xCell) %in% tree.subset$Sample)]
            #valuesOrder <- match(names(values),tree.subset$Sample)
            #values<-values[valuesOrder]
            #tree.subset$Value <- values
            #modelString <- paste("tree.subset$Value + ",modelString,sep="")
            ##########
            
            
            modelString <- paste("~ ",modelString, sep="")
            
            #print(OrganString)
            #print(ImmCellType.xCell)
            #print(modelString)
            tcmCat <- c(OrganString,ImmCellType.xCell,modelString)
            tissCellModel <- rbind(tissCellModel,tcmCat)
            
        }}
    
    
    
    #}# Comment here to just produce the models
    
    
    if(!is.na(numInfl)){
        if(numInfl >= 6 && numNonInfl >= 6)
        {
            print(OrganString)
            print(ImmCellType.xCell)
            intermediateTable <- tree.subset
            rownames(intermediateTable) <- NULL
            intermediateTable$Tissue <- OrganString
            intermediateTable$CellType <- ImmCellType.xCell
            if(p>=185)
            {
                intermediateTable$Autolysis <- NA
                intermediateTable <- intermediateTable[,c(1:3,10,4:9)]
            }
            sampleClustersAndCovariates <- rbind(sampleClustersAndCovariates,intermediateTable)
            
            design <- model.matrix(as.formula(modelString))
            print("Performing Differential Expression...")
            rownames(design) <- tree.subset$Sample
            cols <- colnames(design)
            edgeR.DGElist <- DGEList(counts = geneClustComparison.df,group = tree.subset$Infiltration)###
            keep <- rowSums(cpm(edgeR.DGElist) >= 1) >=5
            edgeR.DGElist <- edgeR.DGElist[keep,]
            edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
            voomTfd <- voom(edgeR.DGElist, design, plot=FALSE)
            voom.fit <- lmFit(voomTfd, design = design)
            voom.fit <- eBayes(voom.fit)
            DGE.results_limma <- topTable(voom.fit,coef = ncol(design),number = Inf, adjust.method = "BH", sort.by = "logFC")
            genes <- rownames(DGE.results_limma)
            DGE.results_limma$Gene <- genes
            DGE.results_limma <- DGE.results_limma[,c(7,1,2,3,4,5,6)]
            print("Writing results...")
            write.table(DGE.results_limma, file = paste(DeconvDGEPath,"PureCell-",OrganString,"-",ImmCellType.xCell,"-Consensus","-DGEresults-ALL.txt",sep=""),quote = FALSE,row.names=FALSE,sep="\t")
            #write.table(DGE.results_limma, file = paste(DeconvDGEPath,"PureCell-",OrganString,"-",ImmCellType.xCell,"-Consensus","-DGEresults-Values.txt",sep=""),quote = FALSE,row.names=FALSE,sep="\t")
        }}
    else
    {
        
        #print(numInflString)
        #print(numNonInflString)
        print("Moving on...")
    }
    
}
write.table(sampleClustersAndCovariates,file=paste(pureCellFilePath,"SampleClusterAssignmentsAndCovariateDataUsedForDGE.txt",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
colnames(tissCellModel) <- c("Tissue","CellType","Design Formula")
write.table(tissCellModel,file = paste(pureCellFilePath,"Tissue-CellType-DesignFormula.txt",sep=""),quote = FALSE,row.names=FALSE,sep="\t")

###################################################################################
###################################################################################
###################################################################################
#DGE post-processing/filtering and preparation for IPA
library(stringr)
pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"
DeconvDGEPath <- paste(pureCellFilePath,"PureCellConsensusDGE",sep="")
fileNames <- read.table(paste(DeconvDGEPath,"/DGEFiles.txt",sep=""),header=F, sep=",")[,1]
DGEFiles <- file.path(DeconvDGEPath,fileNames)
fileNames <- read.table("/Users/ManikUppal/Desktop/ValueRegressed/DGEFiles.txt",header=F, sep=",")[,1]
DGEFiles <- file.path("/Users/ManikUppal/Desktop/ValueRegressed",fileNames)

counter <- 0
masterDF <- data.frame(Gene=character(),logFC=numeric(),Tissue=character(),CellType=character())
for(i in 1:length(DGEFiles))
{
    file <- DGEFiles[i]
    DGETable <- read.table(file,header = TRUE,sep = "\t",quote="")
    logfc.inds <- which(abs(DGETable$logFC) >= 2.0)
    DGETable.lfc <- DGETable[logfc.inds,]
    qval.inds <- which(DGETable.lfc$adj.P.Val < 0.01)
    DGETable.filt <- DGETable.lfc[qval.inds,]
    #print(file)
    
    outputFile <- str_replace_all(file, "ALL", "FILT")
    if(nrow(DGETable.filt) >4 )
    {
        #print(nrow(DGETable.filt))
        #print(max(DGETable.filt$logFC))
        #print(min(DGETable.filt$logFC))
        counter <- counter +1
        write.table(DGETable.filt,outputFile,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
        filename.half <- gsub(paste(DeconvDGEPath,"/PureCell-",sep=""),"",file)
        filename.work <- gsub("-Consensus-DGEresults-ALL.txt","",filename.half)
        #print(filename.work)
        #filename.half <- gsub("/Users/ManikUppal/Desktop/ValueRegressed/PureCell-","",file)
        #filename.work <- gsub("-Consensus-DGEresults-Values.txt","",filename.half)
        
        dashInstances <- unlist(gregexpr(pattern='-',filename.work))
        splitIndex <- 0
        if(length(dashInstances)==1)
        {
            splitIndex <- dashInstances-1
        }
        if(length(dashInstances)==2)
        {
            splitIndex <- dashInstances[2]-1
        }
        Tissue <- substr(filename.work,1,splitIndex)
        CellType <- substr(filename.work,splitIndex+2,nchar(filename.work))
        DGETable.filt$Tissue <- Tissue
        DGETable.filt$CellType <- CellType
        masterDF <- rbind(masterDF, DGETable.filt[,c(1,2,6,8,9)])#can take out the 6 (adj. p value)
        
    }
}
#write.table(masterDF,file=paste(pureCellFilePath,"SignificantDGETranscripts.txt",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
xCellSigs <- read.table("/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/Misc. Intermediate Workup/GeneSignatures/xCellSigs.txt", header=FALSE)
xCellSigs <- as.character(xCellSigs[,1])
CibersortSigs <- read.table("/Users/ManikUppal/Downloads/LM22.txt",header=TRUE,sep="\t")
CibersortSigGenes <- as.character(CibersortSigs[,1])
geneIndsToRemove <- which(masterDF$Gene %in% c(xCellSigs,CibersortSigGenes))
finalGeneDF <- masterDF[-geneIndsToRemove,]
#write.table(finalGeneDF,file=paste(pureCellFilePath,"FilteredDGETranscripts.txt",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
###################################################################################
###################################################################################
###################################################################################
#IPA Results Formatting
#1) Canonical Pathways
#2) Gene Ontology Analysis
#3) Upstream Regulator Analysis
library(stringr)

pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"
IPAFilePath <- paste(pureCellFilePath,"IPAResults/",sep="")

suffix <- c("-Pathways.txt", "-UpstreamRegulators.txt","-GeneOntology.txt")
fileList <- c("CanonicalPathwayFiles.txt","UpstreamRegulatorFiles.txt","GeneOntologyFiles.txt")

for(i in 1:length(suffix))
{
    fileNames <- read.table(paste(IPAFilePath,fileList[i],sep=""),header=F,sep=",")[,1]
    IPAFiles <- file.path(IPAFilePath,fileNames)
    masterIPADF <- data.frame(Sample=character(),Cluster=numeric(),Infiltration=character(),Autolysis=numeric(),Site=character(),SEX=numeric(),AGE=numeric(),DTHHRDY=numeric(),Tissue=character(),CellType=character())
    
    for(j in 1:length(IPAFiles))
    {
        file <- IPAFiles[j]
        filename.half <- gsub(paste(IPAFilePath,"/",sep=""),"",file)
        filename.work <- gsub(suffix[i],"",filename.half)
        dashInstances <- unlist(gregexpr(pattern='-',filename.work))
        splitIndex <- 0
        if(length(dashInstances)==1)
        {
            splitIndex <- dashInstances-1
        }
        if(length(dashInstances)==2)
        {
            splitIndex <- dashInstances[2]-1
        }
        Tissue <- substr(filename.work,1,splitIndex)
        CellType <- substr(filename.work,splitIndex+2,nchar(filename.work))
        IPATable <- read.table(file,header = TRUE,fill = TRUE,skip=2 ,sep = "\t",quote="")
        IPATable$Tissue <- Tissue
        IPATable$CellType <- CellType
        
        if(j>1){
            if(ncol(IPATable) > ncol(masterIPADF))
            {
                missingCols <- colnames(IPATable)[(!(colnames(IPATable) %in% colnames(masterIPADF)))]
                initialColLength <- length(colnames(masterIPADF))
                for(p in 1:length(missingCols))
                {
                    col <- missingCols[p]
                    masterIPADF$New <- NA
                    colnames(masterIPADF)[initialColLength+p] <- col
                }
                seIdx <- match(colnames(IPATable), colnames(masterIPADF))
                masterIPADF <- masterIPADF[,seIdx]
            }
            if(ncol(IPATable) < ncol(masterIPADF))
            {
                missingCols <- colnames(masterIPADF)[(!(colnames(masterIPADF) %in% colnames(IPATable)))]
                initialColLength <- length(colnames(IPATable))
                for(p in 1:length(missingCols))
                {
                    col <- missingCols[p]
                    IPATable$New <- NA
                    colnames(IPATable)[initialColLength+p] <- col
                }
                seIdx <- match(colnames(masterIPADF), colnames(IPATable))
                IPATable <- IPATable[,seIdx]
            }
        }
        
        #print(Tissue)
        #print(CellType)
        #print(ncol(IPATable))
        masterIPADF <- rbind(masterIPADF, IPATable)##We will need to see how to transform this
    }
    write.table(masterIPADF,paste(pureCellFilePath,"SupplementaryTable-IPA",suffix[i],sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
    
}
###################################################################################
###################################################################################
###################################################################################
#Generate the heatmaps for the graphical methods
pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"
DeconvDGEPath <- paste(pureCellFilePath,"PureCellConsensusDGE",sep="")
k_mat <- read.table(paste(pureCellFilePath,"TissCellPairs-k_opt-Updated.txt",sep=""),header=TRUE,stringsAsFactors=FALSE,sep="\t")
k_opt.xCell <- k_mat$`k.xCell`
k_opt.Cibersort <- k_mat$`k.Cibersort.Abs`
k_opt.CibersortRel <- k_mat$`k.Cibersort.Rel`
OrganString <- "Breast - Mammary Tissue"
ImmCellType.Cibersort <- "MastSum"
ImmCellType.xCell <- "Mast cells"
ImmCellType.CibersortRel <-ImmCellType.Cibersort
p <- intersect(which(k_mat$Tissue == OrganString),which(k_mat$xCell == ImmCellType.xCell))


consensusAssignments <- read.table(paste(pureCellFilePath,"HotColdMatrixModels/PureCell-",OrganString,"-",ImmCellType.xCell,"-Consensus","-HotColdMatrix.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")
consensusAssignments <- consensusAssignments[order(consensusAssignments$Cluster),]
hotColdTransitionRow <- table(consensusAssignments$Cluster)[1]
samples.subset <- as.vector(consensusAssignments$Sample)
consensusColdSamps <-samples.subset[1:hotColdTransitionRow]
consensusHotSamps <-samples.subset[(hotColdTransitionRow+1):length(samples.subset)]

k.xCell <- k_opt.xCell[p]
ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.xCell,"-xCell",sep="")
consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.xCell,".consensusClass.csv",sep="")
optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
colnames(optClusters) <- c("Sample","Cluster")
rownames(optClusters) <- optClusters$Sample

coldCluster <- unique(optClusters[which(optClusters$Sample %in% consensusColdSamps),2])
coldSamps <- optClusters[which(optClusters$Cluster == coldCluster),1]
coldSampsNotConsensus <- coldSamps[-which(coldSamps %in% consensusColdSamps)]
finalColdSampOrder <- c(consensusColdSamps,coldSampsNotConsensus)
hotCluster <- unique(optClusters[which(optClusters$Sample %in% consensusHotSamps),2])
hotSamps <- optClusters[which(optClusters$Cluster == hotCluster),1]
hotSampsNotConsensus <- hotSamps[-which(hotSamps %in% consensusHotSamps)]
finalHotSampOrder <- c(hotSampsNotConsensus, consensusHotSamps)

SampleOrder <- c(finalColdSampOrder,finalHotSampOrder)


#load xCell hot and cold results
organ_deconv.xCell <- isolateOrganComparison(deconv_data.xCell,anno_dat,OrganString)
rowIndImm.xCell <- which(rownames(organ_deconv.xCell) %in% ImmCellType.xCell)
ImmCell.organ.xCell <- as.vector(organ_deconv.xCell[rowIndImm.xCell,])
attr(ImmCell.organ.xCell,"names") <- colnames(organ_deconv.xCell)
seIdx <- match(SampleOrder, names(ImmCell.organ.xCell))
ImmCell.organ.xCell <- ImmCell.organ.xCell[seIdx]
triplicate.xCell <- rbind(ImmCell.organ.xCell,ImmCell.organ.xCell,ImmCell.organ.xCell)
q <- pheatmap(triplicate.xCell,annotation_col = optClusters, cluster_rows = FALSE,cluster_cols = FALSE, cex = 0.8, angle_col = 45, main = "Consensus Clustering of Breast Mast cell Expression - xCell", xlab = "Samples")

scale <- 3
tiff("/Users/ManikUppal/Desktop/xCellHeatmap.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(q)
dev.off()

organ_deconv.Cibersort <- isolateOrganComparison(deconv_data.Cibersort,anno_dat,OrganString)
k.Cibersort <- k_opt.Cibersort[p]
rows.Cibersort <- rownames(organ_deconv.Cibersort)
rows.Cibersort <- gsub("\\.", " ",rows.Cibersort)
rowIndImm.Cibersort <- which(rows.Cibersort  %in% ImmCellType.Cibersort)
ImmCell.organ.Cibersort <- as.vector(organ_deconv.Cibersort[rowIndImm.Cibersort,])
attr(ImmCell.organ.Cibersort,"names") <- names(organ_deconv.Cibersort[rowIndImm.Cibersort,])
seIdx <- match(SampleOrder, names(ImmCell.organ.Cibersort))
ImmCell.organ.Cibersort <- ImmCell.organ.Cibersort[seIdx]
triplicate.Cibersort <- rbind(ImmCell.organ.Cibersort,ImmCell.organ.Cibersort,ImmCell.organ.Cibersort)

ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.Cibersort,"-CibersortAbs",sep="")
consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.Cibersort,".consensusClass.csv",sep="")
optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
colnames(optClusters) <- c("Sample","Cluster")
rownames(optClusters) <- optClusters$Sample
q<-pheatmap(triplicate.Cibersort,annotation_col = optClusters, cluster_rows = FALSE,cluster_cols = FALSE, cex = 0.8, angle_col = 45, main = "Consensus Clustering of Breast Mast cell Expression - Cibersort Abs", xlab = "Samples")
scale <- 3
tiff("/Users/ManikUppal/Desktop/CibAbsHeatmap.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(q)
dev.off()

organ_deconv.CibersortRel <- isolateOrganComparison(deconv_data.CibersortRel,anno_dat,OrganString)
k.CibersortRel <- k_opt.CibersortRel[p]
rows.CibersortRel <- rownames(organ_deconv.CibersortRel)
rows.CibersortRel <- gsub("\\.", " ",rows.CibersortRel)
rowIndImm.CibersortRel <- which(rows.CibersortRel  %in% ImmCellType.CibersortRel)
ImmCell.organ.CibersortRel <- as.vector(organ_deconv.CibersortRel[rowIndImm.CibersortRel,])
attr(ImmCell.organ.CibersortRel,"names") <- names(organ_deconv.CibersortRel[rowIndImm.CibersortRel,])
seIdx <- match(SampleOrder, names(ImmCell.organ.CibersortRel))
ImmCell.organ.CibersortRel <- ImmCell.organ.CibersortRel[seIdx]

triplicate.CibersortRel <- rbind(ImmCell.organ.CibersortRel,ImmCell.organ.CibersortRel,ImmCell.organ.CibersortRel)
ConsClustPathName <- paste("PureCell-",OrganString,"-",ImmCellType.CibersortRel,"-CibersortRel",sep="")
consFileName <- paste(pureCellFilePath,"ConsensusClusterResults/",ConsClustPathName,"/",ConsClustPathName,".k=",k.CibersortRel,".consensusClass.csv",sep="")
optClusters <- read.table(consFileName,header=FALSE,sep=",",stringsAsFactors=FALSE)
colnames(optClusters) <- c("Sample","Cluster")
rownames(optClusters) <- optClusters$Sample
subClust1 <- optClusters[which(optClusters$Cluster == 1),]
subClust2 <- optClusters[which(optClusters$Cluster == 3),]
subClust2$Cluster <- 2
subClust3 <- optClusters[which(optClusters$Cluster == 2),]
subClust3$Cluster <- 3
subClust23 <- rbind(subClust2,subClust3)
optClusters <- rbind(subClust1,subClust23)
seIdx <- match(SampleOrder, rownames(optClusters))
optClusters <- optClusters[seIdx,]

q<-pheatmap(triplicate.CibersortRel,annotation_col = optClusters, cluster_rows = FALSE,cluster_cols = FALSE,cex = 0.8, angle_col = 45, main = "Consensus Clustering of Breast Mast cell Expression - Cibersort Rel", xlab = "Samples")
scale <- 3
tiff("/Users/ManikUppal/Desktop/CibRelHeatmap.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(q)
dev.off()

#At this point, use photoshop for figure creation

###################################################################################
###################################################################################
###################################################################################
#DGE figure generation
pureCellFilePath <- "/Users/ManikUppal/Desktop/ElementoLab/GTExPaper/RevisedAnalysis/"
consHotCold <- read.table(paste(pureCellFilePath,"ConsensusHotColdNumbers-New.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
maxSamp <- max(consHotCold[which(!is.na(consHotCold[,3]+consHotCold[,4])),3] +consHotCold[which(!is.na(consHotCold[,3]+consHotCold[,4])),4])
finalGeneDF <- read.table(paste(pureCellFilePath,"FilteredDGETranscripts.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
#targetGenesList <- as.character(read.table("/Users/ManikUppal/Desktop/DGEFilesFilt.txt",header=F, fill=TRUE,sep=",")[,1])

cellTypes.df <- data.frame(ciber=c('T cells CD8','T cells CD4 naive','CD4_memory','Neutrophils','MacrophageSum','Bcellsum','NK_Sum','DendriticSum','MastSum','Myeloid_Sum','T cells follicular helper','T cells regulatory (Tregs)','T cells gamma delta','Monocytes','Eosinophils','Lymph_Sum'), xcell=c('CD8Sum','CD4+ naive T-cells','CD4_memory','Neutrophils','MacrophageSum','Bcellsum','NK cells','DendriticSum','Mast cells','Myeloid_Sum','Th_Sum','Tregs','Tgd cells','Monocytes','Eosinophils','Lymph_Sum'),stringsAsFactors = F)

#cellTissTable <- matrix(,nrow=maxSamp,ncol=length(targetGenesList))
names <- character(0)
LeDataFrame <- data.frame(Infiltration = double(),Sample=character(),Cluster=integer(),TissCell=character())
for(i in 1:nrow(unique(finalGeneDF[,c(4,5)])))
{
    Tissue <- unique(finalGeneDF[,c(4,5)])[i,1]
    CellType <- unique(finalGeneDF[,c(4,5)])[i,2]
    tissCellModelMat <- read.table(paste(pureCellFilePath,"HotColdMatrixModels/PureCell-",Tissue,"-",CellType,"-Consensus-HotColdMatrix.txt",sep=""),header=TRUE,sep="\t")
    
    name <- paste(Tissue, CellType,sep="-")
    names <-c(names,name)
    CellType.ciber <- cellTypes.df[which(cellTypes.df$xcell == CellType),1]
    CellType <- CellType.ciber
    if(CellType == "T cells CD8")
    {
        CellType <- "T.cells.CD8"
    }
    if(CellType == "T cells follicular helper")
    {
        CellType <- "T.cells.follicular.helper"
    }
    cell.type.ind <- which(rownames(deconv_data.Cibersort) == CellType)
    organ_deconv <- isolateOrganComparison(deconv_data.Cibersort, anno_dat, Tissue)
    inds <- which(colnames(organ_deconv) %in% tissCellModelMat$Sample)
    #organ_deconv.cell <- organ_deconv[cell.type.ind,]
    #if(i==26)
    #{
    #organ_deconv.cell[407] <- as.numeric(organ_deconv.cell[407])
    #}
    #organ_deconv.cell.scld <- scale(t(organ_deconv.cell))
    #organ_deconv.rd <- as.vector(organ_deconv.cell.scld[inds,])
    organ_deconv.cell <- as.matrix(organ_deconv[cell.type.ind,])
    organ_deconv.cell.scld <- scale(organ_deconv.cell)
    organ_deconv.rd <- as.vector(organ_deconv.cell.scld[inds,])
    organ_deconv.test <- as.data.frame(organ_deconv.rd)
    organ_deconv.test$Sample <- names(organ_deconv.cell[inds,])
    seIdx <- match(tissCellModelMat$Sample, organ_deconv.test$Sample)
    organ_deconv.test <- organ_deconv.test[seIdx,]
    organ_deconv.test$Cluster <- tissCellModelMat$Infiltration
    organ_deconv.test$TissCell <- name
    organ_deconv.test$Tissue <- Tissue#
    organ_deconv.test$CellType <- CellType#
    #print(paste(Tissue,CellType,"Scaled",mean(organ_deconv.cell.scld),sd(organ_deconv.cell.scld),"Reduced",mean(organ_deconv.rd),sd(organ_deconv.rd),sep=" "))
    ##length(organ_deconv.rd) <- nrow(cellTissTable)
    ##cellTissTable[,i] <- organ_deconv.rd
    LeDataFrame <- rbind(LeDataFrame,organ_deconv.test)
    
}
#rownames(LeDataFrame) <- LeDataFrame$Sample
LeDataFrame <- LeDataFrame[,c(2,1,4,3)]
colnames(LeDataFrame) <- c("Sample","Infiltration","TissCell","Cluster")
LeDataFrame$Cluster <- factor(LeDataFrame$Cluster)
p <- ggplot(LeDataFrame, aes(x=TissCell,y=Infiltration,fill=Cluster)) + geom_boxplot(outlier.shape=NA) + geom_jitter(size=0.05,position=position_jitter(width=.05, height=0)) + geom_point(aes(color=Cluster)) + geom_point(shape = 1,size = .05,colour = "black") + labs(x = "Tissue-Cell Pairs", y = "Z-Score of Infiltration value\n(Cibersort Absolute)") + theme(text = element_text(size=5),axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 15, face = "bold"), panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey"), plot.margin = margin(1,0,1,2,"cm")) + scale_fill_manual(values=c("slateblue1","firebrick")) + scale_color_manual(values=c("blue","red"))

scale <- 3
tiff("/Users/ManikUppal/Desktop/HotColdClusters-New.png",width=900*scale,height=600*scale,res=100*scale,units="px")
print(p)
dev.off()



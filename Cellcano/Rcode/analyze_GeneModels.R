## RSPH Seurat environment
## Try a different naming scheme
workDir = "~/data/GSE129785_scATACseq"
setwd(workDir)

modelDir = file.path(workDir, "ArchR_GeneModels")
resultDir = file.path(workDir, "ArchR_GeneModels_analysis")
dir.create(resultDir, showWarnings=F)

library(Matrix)

## ======== Configuration ========
randomized = F ## if randomly selected 500 genes along with 1000 cells
matrixNorm = F
## ===============================

modelFiles = list.files(modelDir, pattern="*.mtx", full.names=T)
#modelFiles = list.files(modelDir, pattern="*.mtx", full.names=T)[1:3]  ## for testing

allGenes = lapply(seq(modelFiles), function(i) {
    modelFile = modelFiles[i]
    modelGenes = scan(gsub('.mtx', '_genes.tsv', modelFile), what=character())
    modelGenes
})

allCells = lapply(seq(modelFiles), function(i) {
    modelFile = modelFiles[i]
    modelCells = scan(gsub('.mtx', '_cells.tsv', modelFile), what=character())
    modelCells
})

commonGenes = Reduce(intersect, allGenes) ## get common genes
commonCells = Reduce(intersect, allCells) ## get common genes

if (randomized) {
set.seed(2021)
sampledGenes = sample(commonGenes, size=500)
set.seed(2021)
sampledCells = sample(commonCells, size=1000)

modelGeneList = lapply(seq(modelFiles), function(i) {
    modelFile = modelFiles[i]
    model = readMM(modelFile)
    modelGenes = scan(gsub('.mtx', '_genes.tsv', modelFile), what=character())
    modelCells = scan(gsub('.mtx', '_cells.tsv', modelFile), what=character())
    rownames(model) = modelGenes
    colnames(model) = modelCells

    ## extract model from sample
    #sampleCells = colnames(model)[grepl(sample, colnames(model))]
    sampleMat = model[sampledGenes, sampledCells]  ## extract specific gene
    as.vector(sampleMat)
})
names(modelGeneList) = gsub('.mtx', '', basename(modelFiles))
saveRDS(modelGeneList, file=file.path(resultDir, "sampledGeneList.RDS"))

## generate correlation matrix
corMat = matrix(0, length(modelGeneList), length(modelGeneList))
rownames(corMat) = names(modelGeneList)
colnames(corMat) = names(modelGeneList) 

for (modelGene1 in names(modelGeneList)) {
    for (modelGene2 in names(modelGeneList)) {
        corMat[modelGene1, modelGene2] = cor(modelGeneList[[modelGene1]], modelGeneList[[modelGene2]], method="spearman")
    }
}
saveRDS(corMat, file.path(resultDir, "SpearmanCor_sampledGeneModels.RDS"))

## draw out
library(reshape2)
library(ggplot2)

meltedDf = melt(corMat)
colnames(meltedDf) = c("GeneModel1", "GeneModel2", "Corr")
ggplot(meltedDf, aes(GeneModel1, GeneModel2, fill=Corr)) + 
    geom_tile(color="black") +
    coord_fixed() + ## square tile
    theme(axis.text.x=element_text(size=6, angle=90, vjust=1, hjust=1),
          axis.text.y=element_text(size=6)) + 
    scale_fill_gradient(low="white", high="red")
 
ggsave(file.path(resultDir, "SpearmanCor_sampledGeneModels.png"))
}

if (matrixNorm) {
## may not be needed
modelGeneList = lapply(seq(modelFiles), function(i) {
    modelFile = modelFiles[i]
    model = readMM(modelFile)
    modelGenes = scan(gsub('.mtx', '_genes.tsv', modelFile), what=character())
    modelCells = scan(gsub('.mtx', '_cells.tsv', modelFile), what=character())
    rownames(model) = modelGenes
    colnames(model) = modelCells

    ## extract model from sample
    #sampleCells = colnames(model)[grepl(sample, colnames(model))]
    mat = model[commonGenes, commonCells]  ## extract specific gene
    mat
})
names(modelGeneList) = gsub('.mtx', '', basename(modelFiles))
saveRDS(modelGeneList, file=file.path(resultDir, "allGeneList.RDS"))
}


modelGeneList = readRDS(file.path(resultDir, "allGeneList.RDS"))
set.seed(0)
sampledCells = sample(commonCells, size=100)
## perform PCA on one selected cell
cellPCAList = lapply(seq(sampledCells), function(i) {
    cellName = commonCells[i]
    cat(i, cellName, '\n')
    cellGeneList = lapply(seq(modelGeneList), function(j) {
        cellGene = modelGeneList[[j]]
        cellGene[, cellName]
    })
    ## turn into dataframe and conduct PCA analysis
    cellGeneDF = Reduce(cbind, cellGeneList)
    colnames(cellGeneDF) = 1:ncol(cellGeneDF)
    pca_res = prcomp(cellGeneDF, center=T, scale=T)
    pca_summary = summary(pca_res)
    return(pca_summary)
})
names(cellPCAList) = names(commonCells)
saveRDS(cellPCAList, file=file.path(resultDir, "sampled_cellwisePCA.RDS"))

cellPCAList = readRDS(file.path(resultDir, "sampled_cellwisePCA.RDS"))
cellPC1toPC5Variance = lapply(seq(cellPCAList), function(i) {
    cellPCA = cellPCAList[[i]]
    cellPCA$importance[3, 1:5]
})
cellPC1toPC5VarianceDf = Reduce(rbind, cellPC1toPC5Variance)
rownames(cellPC1toPC5VarianceDf) = NULL

library(reshape2)
library(ggplot2)
meltedDf = melt(cellPC1toPC5VarianceDf)
g = ggplot(meltedDf, aes(x=Var2, y=value)) + 
    geom_boxplot() +
    ylim(c(0.65, 1)) +
    labs(x = "PC", y = "Cumulative Variances Explained")
ggsave(file.path(resultDir, "sampled_cellwisePCA_boxplot.png"))


#cellName = commonCells[1]
#cellGeneList = lapply(seq(modelGeneList), function(i) {
#    cellGene = modelGeneList[[i]]
#    cellGene[, cellName]
#})
#saveRDS(cellGeneList, file=file.path(resultDir, "oneCellGeneList.RDS"))

#cellGeneList = readRDS(file.path(resultDir, "oneCellGeneList.RDS"))


## write as h5df, jobs kept getting killed
#library(rhdf5)
#library(HDF5Array)
#h5fileDir = file.path(resultDir, "allGeneList.h5")
#h5createFile(h5fileDir)
#lapply(seq(modelGeneList), function(i) {
#    mat = modelGeneList[[i]]
#    groupName = names(modelGeneList)[i]
#    writeHDF5Array(mat, filepath=h5fileDir, name=groupName)
#})
#h5closeAll()

## apply lognorm + scale, TODO: memory issue!
#scaledModels = lapply(seq(modelGeneList), function(i) {
#    X = modelGeneList[[i]]
#    k = colSums(X)
#    k = k/median(k)
#    Xnorm = sweep(X, 2, k, FUN="/")
#    Xlognorm = log(1+Xnorm)
#    Xscaled = scale(t(Xlognorm))
#    Xscaled[is.nan(Xscaled)] = 0
#    Xscaled
#})
#saveRDS(scaledModels, file=file.path(resultDir, "scaledGeneList.RDS"))

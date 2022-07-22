work_dir = "~/data/GSE129785_scATACseq"
setwd(work_dir)

model_dir = file.path(work_dir, "ArchR_GeneModels")
dir.create(model_dir, showWarnings=F)

library(ArchR)
set.seed(1)
addArchRThreads(threads=4) 
addArchRGenome("hg19")  ## according to data, it was sequenced in hg19
genes <- getGenes()
mcols(genes)$name <- mcols(genes)$symbol

proj = loadArchRProject(path = "PBMC")

## model promoter regions
dropStrand <- function(gr){
	strand(gr) <- "*"
	gr
}

## extract matrix info
save_genescore <- function(matObj, modelName) {
    mat = assays(matObj)[[1]]
    mat_barcodes = sapply(strsplit(colnames(mat), split='#'), '[', 2)
    mat_rep = substr(colnames(mat), 12, 20)
    colnames(mat) = paste(mat_rep, mat_barcodes, sep='#')
    rowMeta = rowData(matObj)
    rowIdx = which(!duplicated(rowMeta$name))
    mat = mat[rowIdx, ]
    rownames(mat) = rowMeta[rowIdx, ]$name
    writeMM(mat, file.path(model_dir, paste0(modelName, '.mtx')))
    write(rownames(mat), file.path(model_dir, paste0(modelName, '_genes.tsv')))
    write(colnames(mat), file.path(model_dir, paste0(modelName, '_cells.tsv')))
}

## add 54 gene score models

## ====== Promoter regions
i = 1

cat(i, "Adding Model-Promoter-1...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features=dropStrand(resize(resize(genes, 1, "start"), 1000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-1")
i = i + 1

cat(i, "Adding Model-Promoter-2...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 2000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-2")
i = i + 1

cat(i, "Adding Model-Promoter-5...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 5000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-5")
i = i + 1

cat(i, "Adding Model-Promoter-10...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 10000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-10")
i = i + 1

cat(i, "Adding Model-Promoter-25...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 25000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-25")
i = i + 1

cat(i, "Adding Model-Promoter-50...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 50000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-50")
i = i + 1

cat(i, "Adding Model-Promoter-100...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(resize(resize(genes, 1, "start"), 100000, "center")),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-Promoter-100")
i = i + 1

## ===== Gene body extension
cat(i, "Adding Model-GeneBody-0_0...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 0, 0)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-0_0")
i = i + 1

cat(i, "Adding Model-GeneBody-1000_0...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 1000, 0)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-1000_0")
i = i + 1

cat(i, "Adding Model-GeneBody-2000_0...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 2000, 0)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-2000_0")
i = i + 1

cat(i, "Adding Model-GeneBody-5000_0...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 5000, 0)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-5000_0")
i = i + 1

cat(i, "Adding Model-GeneBody-10000_0...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 10000, 0)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-10000_0")
i = i + 1

cat(i, "Adding Model-GeneBody-1000_1000...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 1000, 1000)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-1000_1000")
i = i + 1

cat(i, "Adding Model-GeneBody-2000_2000...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 2000, 2000)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-2000_2000")
i = i + 1

cat(i, "Adding Model-GeneBody-5000_5000...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 5000, 5000)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix = "FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-5000_5000")
i = i + 1

cat(i, "Adding Model-GeneBody-10000_10000...\n")
addFeatureMatrix(input=proj, matrixName="FeatureMatrix",
                 features = dropStrand(extendGR(genes, 10000, 10000)),
                 force=T)
matObj = getMatrixFromProject(proj, useMatrix="FeatureMatrix")
save_genescore(matObj, modelName="Model-GeneBody-10000_10000")
i = i + 1

### ========= Gene Models From TSS
j = 1  ## for counting gene model
cat(i, paste0("Adding GeneModel-Constant-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
        geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 5000),
	extendDownstream = c(1000, 5000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
) 
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-Constant-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-Constant-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 10000),
	extendDownstream = c(1000, 10000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5,
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-Constant-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-Constant-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 25000),
	extendDownstream = c(1000, 25000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-Constant-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-Constant-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "1",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-Constant-", j))
i = i + 1
j = j + 1

### ========= different Gene Models From TSS
j = 1  ## for counting gene model

cat(i, paste0("Adding GeneModel-TSS-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-Exponential-", j))
i = i + 1
j = j + 1


### =========  Gene Models From TSS No Boundary
j = 1  ## for counting gene model

cat(i, paste0("Adding GeneModel-TSS-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-TSS-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = TRUE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-TSS-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

### =========  Custom Gene Models From GB
j = 1  ## for counting gene model

cat(i, paste0("Adding GeneModel-GB-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-", j))
i = i + 1
j = j + 1

### =========  Custom Gene Models From TSS No Boundary
j = 1  ## for counting gene model
cat(i, paste0("Adding GeneModel-GB-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1


cat(i, paste0("Adding GeneModel-GB-NoBoundary-Exponential-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/100000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 0, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = FALSE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-NoBoundary-Exponential-", j))
i = i + 1
j = j + 1

### =========  Custom Gene Models From GB extend
j = 1  ## for counting gene model
cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/10000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
	genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/25000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 0, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 1000, #New Param
	geneDownstream = 1000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 2000, #New Param
	geneDownstream = 2000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))
i = i + 1
j = j + 1

cat(i, paste0("Adding GeneModel-GB-Exponential-Extend-", j), '...\n')
addGeneScoreMatrix(input=proj,
        genes = genes,
	geneModel = "exp(-abs(x)/5000) + exp(-1)",
	matrixName = "GeneScoreMatrix",
	extendUpstream = c(1000, 100000),
	extendDownstream = c(1000, 100000),
	geneUpstream = 5000, #New Param
	geneDownstream = 5000, #New Param
	useGeneBoundaries = TRUE,
	useTSS = FALSE, #New Param
	tileSize = 500,
	ceiling = 4,
	geneScaleFactor = 5, #New Param
        force = T
)
matObj = getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
save_genescore(matObj, modelName=paste0("GeneModel-GB-Exponential-Extend-", j))


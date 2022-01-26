work_dir = "/projects/compbio/users/wma36/data/10X_Multiome_PBMC10K"
setwd(work_dir)

library(ArchR)
set.seed(1)
addArchRThreads(threads = 8) 
addArchRGenome("hg19")  ## according to data, it was sequenced in hg38, but i lifted over to hg19


## load proj
#proj = loadArchRProject(path="10XMultiome_PBMC")

## need to be bgzip files not gzip files
ArrowFiles <- createArrowFiles(
                    inputFiles = file.path(work_dir, "pbmc_granulocyte_sorted_10k_atac_fragments.hg19.sorted.tsv.gz"),
                    sampleNames = "10XMultiome_PBMC",
                    minTSS = 8, #Dont set this too high because you can always increase later
                    ## paper set it to be 8, but it was based on 25M cells, here we have 2M
                    minFrags = 1000, 
                    addTileMat = TRUE,
                    addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
                    input = ArrowFiles,
                    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
                    LSIMethod = 1
)

proj <- ArchRProject(
          ArrowFiles = ArrowFiles, 
          outputDirectory = "10XMultiome_PBMC",
          copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
saveArchRProject(ArchRProj = proj, outputDirectory = "10XMultiome_PBMC", load = FALSE)

## get ColData
proj$UMAP1 = proj@embeddings[["UMAP"]][["df"]][, 1]
proj$UMAP2 = proj@embeddings[["UMAP"]][["df"]][, 2]
g = plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggsave("ArchR_umap.png", g)

PBMC_cells = read.table(file.path(work_dir, "sample_metadata.csv"), header=T, sep=',')

## add doublet scores
proj_df = getCellColData(proj)
proj_df$barcode = sapply(strsplit(rownames(proj_df), '#'), '[', 2)
proj_df = merge(proj_df, PBMC_cells, by.x=c("barcode"), by.y=c("barcode"), all.x=T)
proj_df = proj_df[!is.na(proj_df$celltype), ] ## 10,022 cells
#write.table(as.data.frame(proj_df), "ArchR_cells_df.tsv", row.names=F, sep='\t', quote=F)

library(ggplot2)
ggplot(as.data.frame(proj_df), aes(x=UMAP1, y=UMAP2, color=celltype)) + geom_point(alpha=0.5) + theme(legend.position="bottom")
ggsave("ArchR_umap_celltype.png")

ggplot(as.data.frame(proj_df), aes(x=UMAP1, y=UMAP2, color=DoubletScore)) + geom_point(alpha=0.5) + theme(legend.position="bottom")
ggsave("ArchR_umap_doublet.png")

## add Gene Scores
#proj <- addImputeWeights(proj)  ## smooth dropout in gene scores, I will let it be sparse by now
filtered_proj = proj[paste(proj_df$Sample, proj_df$barcode, sep='#')]
genescore_mat = getMatrixFromProject(filtered_proj, useMatrix = "GeneScoreMatrix")
writeMM(assays(genescore_mat)$GeneScoreMatrix, "ArchR_genescore.mtx")
write(rowData(genescore_mat)$name, "ArchR_genescore_genes.tsv")
genescore_mat_colData = colData(genescore_mat)
genescore_mat_colData$barcodes = sapply(strsplit(rownames(genescore_mat_colData), '#'), '[', 2)
write(genescore_mat_colData$barcodes, "ArchR_genescore_barcodes.tsv")

## extract tile matrix
tile_mat = getMatrixFromProject(filtered_proj, useMatrix="TileMatrix", binarize=T)
writeMM(assays(tile_mat)$TileMatrix, "10XMultiome_tile.mtx")
tile_mat_colData = colData(tile_mat)
tile_mat_colData$barcodes = sapply(strsplit(rownames(tile_mat_colData), '#'), '[', 2)
write(tile_mat_colData$barcodes,  "10XMultiome_tile_barcodes.tsv")


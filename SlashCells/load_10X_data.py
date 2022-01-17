import os
import anndata
import numpy as np
import pandas as pd

## === 10X Multiome PBMC
def curate_10XPBMC_celltypes(adata, celltype_gran=1, celltype_label="cell.type"):
    '''Curate PBMC cell types into common cell types between immune PBMC, AML PBMC and 10X multiome PBMC

    @celltype_gran: default 1 because the given cell.type are sub-cell types; if 0, then curate to major cell types
    '''
    celltype_categories = {'B cells': ['memory B cells', 'naive B cells'],
            'CD4 T cells': ['memory CD4 T cells', 'naive CD4 T cells'],
            'CD8 T cells': ['effector CD8 T cells', 'naive CD8 T cells'],
            'NK cells': ['CD56 (bright) NK cells', 'CD56 (dim) NK cells'],
            'CD14+ Monocytes': ['classical monocytes'],
            'FCGR3A+ Monocytes': ['non-classical monocytes'],
            'Dendritic cells': ['myeloid DC', 'plasmacytoid DC']}
    ## extract out major cell types
    major_celltypes = []
    for celltype in adata.obs[celltype_label].tolist():
        flag = False
        for major_celltype, sub_celltypes in celltype_categories.items():
            if celltype in sub_celltypes:
                major_celltypes.append(major_celltype)

                ## if found
                flag = True
                break
        if False == flag:
            major_celltypes.append(np.nan)
    adata.obs['majortypes'] = major_celltypes

    if 0 == celltype_gran:
        adata.obs.rename(columns={celltype_label: "subtypes", 
                                    "majortypes": celltype_label}, 
                        inplace=True)
    return adata

## ==== 10X PBMC Multiome, sorted, 10K
def add_cell_metadata(data_dir, adata):
    ## load cell metadata
    cell_metadata = pd.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/sample_metadata.csv',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['barcode']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="barcode", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'celltype': 'cell.type'}, inplace=True)
    return adata


def load_10XPBMC_ATACbins(data_dir):
    '''Loading 10X PBMC ATAC bins called by ArchR
    '''
    matrice_dir = data_dir+os.sep+'10X_Multiome_PBMC10K/matrices'
    adata = anndata.read_mtx(matrice_dir+os.sep+'10XMultiome_tile.mtx').T
    cells = pd.read_csv(matrice_dir+os.sep+'10XMultiome_tile_barcodes.tsv',
            header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None
    ## add cell metadata
    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_genescore(data_dir, model_name="GeneModel-GB-Exponential-Extend-18"):
    '''Loading 10X PBMC gene score called by ArchR
    '''
    matrice_dir = data_dir+os.sep+'10X_Multiome_PBMC10K/ArchR_GeneModels'
    ## load adata
    adata = anndata.read_mtx(matrice_dir+os.sep+model_name+'.mtx.gz').T
    ## load cells and genes
    genes = pd.read_csv(matrice_dir+os.sep+model_name+'_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(matrice_dir+os.sep+model_name+'_cells.tsv', header=None, sep='\t')
    adata.obs["barcode"] = [x.replace('10XMultiome_PBMC#', '') for x in cells[0].values.tolist()]
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    ## add cell metadata
    adata = add_cell_metadata(data_dir, adata)
    return adata
 
def load_10XPBMC_ge(data_dir):
    '''Loading 10X PBMC gene expression
    '''
    ## load adata
    adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/filtered_feature_bc_matrix/RNA.mtx').T
    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/filtered_feature_bc_matrix/RNA_features.tsv",
            header=None, sep='\t')

    ## filter gene symbols
    features = pd.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/filtered_feature_bc_matrix/features.tsv',
            header=None, sep='\t')
    gene_df = genes.merge(features, left_on=0, right_on=0, how="left")
    dedup_gene_df = gene_df.drop_duplicates(subset=[1], keep='first', inplace=False)
    adata = adata[:, list(dedup_gene_df.index)]
    adata.var["genes"] = dedup_gene_df[1].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    cells = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/filtered_feature_bc_matrix/barcodes.tsv",
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_ge_SC2P(data_dir):
    '''Loading 10X PBMC SC2P transformed result
    '''

    ## load adata
    adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/SC2P.mtx').T
    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/SC2P_genes.tsv",
            header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    cells = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/SC2P_barcodes.tsv",
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_ArchRgenescore(data_dir):
    '''Loading 10X PBMC ArchR gene score
    '''
    matrice_dir = data_dir+os.sep+'10X_Multiome_PBMC10K/matrices'
    ## load adata
    adata = anndata.read_mtx(matrice_dir+os.sep+'ArchR_genescore.mtx').T
    ## load cells and genes
    genes = pd.read_csv(matrice_dir+os.sep+'ArchR_genescore_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(matrice_dir+os.sep+'ArchR_genescore_barcodes.tsv', header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_cisTopic(data_dir):
    '''Loading 10X PBMC cisTopic transformation with 8 topics
    '''
    adata = anndata.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/cisTopic_result.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_SnapATAC(data_dir, option="raw"):
    '''Loading 10X PBMC SnapATAC

    @option: 
        - raw: bin-by-cell matrix, 5KB window
        - PC: SVD first 16 PCs on bin-by-cell matrix
    '''

    ## load adata
    if option == "PC":
        adata = anndata.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/SnapATAC_16PCs_mat.csv')
        ## load genes and cells
        adata.obs['barcode'] = adata.obs.index
    if option == "raw":
        adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/10X_PBMC_5KB.mtx')
        cells = pd.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/10X_PBMC_5KB_barcodes.tsv',
                header=None, sep='\t')
        adata.obs['barcode'] = cells[0].values
        adata.obs_names = cells[0]

    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs_names_make_unique(join="-")

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_peak(data_dir):
    '''Loading 10X PBMC peak-by-cell count

    '''
    ## load adata
    adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/filtered_feature_bc_matrix/ATAC.mtx').T
    ## load cells and peaks
    peaks = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/filtered_feature_bc_matrix/ATAC_features.tsv",
            header=None, sep='\t')
    adata.var["peaks"] = peaks[0].values
    adata.var_names = adata.var["peaks"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    cells = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/filtered_feature_bc_matrix/barcodes.tsv",
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_promoter(data_dir, binarize=False):
    '''Loading 10X PBMC promoter counts

    @ binarize: whether to binarize the promoter counts
    '''
    ## load data as anndata
    if binarize:
        adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/promoter_binarized_counts.mtx').T
    else:
        adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/promoter_counts.mtx').T
    ## load genes
    genes = pd.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/promoter_genes.tsv',
            header=None, sep='\t')
    adata.var['gene_symbols'] = genes[0].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    ## load cells
    cells = pd.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/promoter_barcodes.tsv',
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_ga(data_dir):
    '''Loading cicero gene activities
    '''
    adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/cicero_gene_activities.mtx').T
    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/cicero_genes.tsv",
            header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    cells = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/cicero_barcodes.tsv",
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_Cusanovich(data_dir):
    ## load 10X PBMC Cusanovich, 20 SVD components
    adata = anndata.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/cusanovich_LSI_matrix.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")
    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_Seurat_imputed(data_dir):
    ## load 10X PBMC Seurat imputed scATAC-seq data
    adata = anndata.read_csv(data_dir+os.sep+'10X_Multiome_PBMC10K/Seurat_imputed.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    barcodes = [x.replace('.', '-') for x in adata.obs.index]  ## replace . with -
    adata.obs.index = barcodes
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")
    adata = add_cell_metadata(data_dir, adata)
    return adata

def load_10XPBMC_Seurat_NNweights(data_dir, opt="WKNN"):
    '''Loading 10X PBMC WKNN/WSNN

    @opt: WKNN/WSNN
    '''
    if opt == "WKNN":
        adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/Seurat_WKNN.mtx').T
    if opt == "WSNN":
        adata = anndata.read_mtx(data_dir+os.sep+'10X_Multiome_PBMC10K/Seurat_WSNN.mtx').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
 
    cells = pd.read_csv(data_dir+os.sep+"10X_Multiome_PBMC10K/Seurat_WKNN_barcodes.tsv",
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None
    adata = add_cell_metadata(data_dir, adata)
    return adata

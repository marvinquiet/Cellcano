import os
import anndata
import numpy as np
import pandas as pd

## === Greenleaf immune study PBMC GSE129785
def curate_immunePBMC_celltypes(adata, celltype_gran=1, celltype_label="cell.type"):
    '''Curate PBMC cell types into common cell types between immune PBMC, AML PBMC and 10X multiome PBMC

    @celltype_gran: default 1 because the given cell.type are sub-cell types; if 0, then curate to major cell types
    '''
    celltype_categories = {'B cells': ['Memory B', 'Naive B', 'Plasma B'],
            'CD4 T cells': ['N CD4 T1', 'N CD4 T2', 'M CD4 T'],
            'CD8 T cells': ['CM CD8 T', 'EM CD8 T', 'N CD8 T1', 'N CD8 T2', 'N CD8 T3'],
            'NK cells': ['Imm NK', 'Mat NK1', 'Mat NK2'],
            'CD14+ Monocytes': ['Mono2'],
            'FCGR3A+ Monocytes': ['Mono1'],
            'Dendritic cells': ['cDC', 'pDC']}
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

def curate_AMLPBMC_celltypes(adata, celltype_gran=1, celltype_label="cell.type"):
    '''Curate PBMC cell types into common cell types between immune PBMC, AML PBMC and 10X multiome PBMC

    @celltype_gran: default 1 because the given cell.type are sub-cell types; if 0, then curate to major cell types
    '''
    celltype_categories = {'B cells': ['B', 'Plasma'],
            'CD4 T cells': ['CD4.N1', 'CD4.N2', 'CD4.M'],
            'CD8 T cells': ['CD8.N', 'CD8.EM', 'CD8.CM'],
            'NK cells': ['NK'],
            'CD14+ Monocytes': ['CD14.Mono.1', 'CD14.Mono.2'],
            'FCGR3A+ Monocytes': ['Unk'],
            'Dendritic cells': ['cDC', 'pDC']}
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



def add_immunePBMC_metadata(data_dir, adata):
    '''Loading cell metadata for immune PBMC 
    '''
    metadata = pd.read_csv(data_dir+os.sep+'GSE129785_PBMC_scATACseq/GSE129785_scATAC-Hematopoiesis-PBMC.cell_barcodes.txt',
            header=0, sep='\t')
    metadata["ID"] = metadata["Group_Barcode"]+'#'+metadata["Barcodes"]

    common_barcodes = set(metadata["ID"]).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(metadata, left_on="barcode", right_on="ID", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'celltype': 'cell.type'}, inplace=True)
    return adata

def load_immunePBMC_genescore(data_dir, ind=None, model_name="GeneModel-GB-Exponential-Extend-18"):
    '''Load immune PBMC gene score matrix derived from ArchR
    '''
    matrice_dir = data_dir+os.sep+'GSE129785_PBMC_scATACseq/ArchR_GeneModels'
    ## load adata
    adata = anndata.read_mtx(matrice_dir+os.sep+model_name+'.mtx.gz').T
    ## load cells and genes
    genes = pd.read_csv(matrice_dir+os.sep+model_name+'_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(matrice_dir+os.sep+model_name+'_cells.tsv', header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_immunePBMC_metadata(data_dir, adata)
    ## filter by individual
    if ind != "all" and ind is not None:
        ind_cells = adata.obs[adata.obs["Group"].isin(ind.split('#'))].index
        adata = adata[ind_cells]
 
    return adata

def load_immunePBMC_ATACbins(data_dir, ind=None):
    '''Load immune PBMC ATAC-seq 500bp bins matrix derived from ArchR

    @ind: replicates separated by # separator (same as how the ArchR did for barcode), can be one of PBMC_Rep1-4
    '''
    matrice_dir = data_dir+os.sep+'GSE129785_PBMC_scATACseq/matrices'

    reps = ["PBMC_Rep1", "PBMC_Rep2", "PBMC_Rep3", "PBMC_Rep4"]  ## 4 individual
    if ind == "all" or ind is None:
        inds = reps
    else:
        inds = ind.split('#')

    adata_list = []
    for rep in inds:
        rep_adata = anndata.read_mtx(matrice_dir+os.sep+rep+'_tile.mtx').T
        rep_cells = pd.read_csv(matrice_dir+os.sep+rep+'_tile_barcodes.tsv', 
                header=None, sep='\t')
        rep_adata.obs["barcode"] = rep_cells[0].values
        rep_adata.obs_names = rep_adata.obs['barcode']
        rep_adata.obs_names_make_unique(join="-")
        rep_adata.obs.index.name = None
        adata_list.append(rep_adata)

    from anndata import AnnData
    adata = AnnData.concatenate(*adata_list, index_unique=None, batch_categories=inds)
    ## add cell metadata
    adata = add_immunePBMC_metadata(data_dir, adata)
    return adata

def load_immunePBMC_ciceroGA(data_dir, ind=None):
    '''Load immune PBMC gene activities derived from cicero

    @ind: replicates separated by # separator (same as how the ArchR did for barcode), can be one of PBMC_Rep1-4
    '''
    matrice_dir = data_dir+os.sep+'GSE129785_PBMC_scATACseq/matrices'
    ## load adata
    adata = anndata.read_mtx(matrice_dir+os.sep+'cicero_gene_activities.mtx').T
    ## load cells and genes
    genes = pd.read_csv(matrice_dir+os.sep+'cicero_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(matrice_dir+os.sep+'cicero_barcodes.tsv', header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_immunePBMC_metadata(data_dir, adata)
    ## filter by individual
    if ind != "all" and ind is not None:
        ind_cells = adata.obs[adata.obs["Group"].isin(ind.split('#'))].index
        adata = adata[ind_cells]
    return adata

def load_immunePBMC_SnapATAC(data_dir, harmony=False, ind=None):
    '''Load immune PBMC dataset after SnapATAC dimension projection

    The projection uses PBMC_Rep1 as reference and do diffusion map, then project PBMC_Rep2 towards it
    @harmony: whether to remove batch effect in dimension reduction space
    @ind: PBMC_Rep1 and PBMC_Rep2 (we only have these two)
    '''

    matrice_dir = data_dir+os.sep+'GSE129785_PBMC_scATACseq/matrices'
    ## load adata
    harmony_str = "after" if harmony else "before"
    adata = anndata.read_mtx(matrice_dir+os.sep+'SnapATAC_'+harmony_str+'_harmony.mtx')
    ## load barcodes
    cells = pd.read_csv(matrice_dir+os.sep+'SnapATAC_'+harmony_str+'_harmony_barcodes.tsv',
            header=None, sep='\t')
    adata.obs["barcode"] = cells[0].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_immunePBMC_metadata(data_dir, adata)
    if ind != "all" or ind is not None:
        ind_cells = adata.obs[adata.obs["Group"].isin(ind.split('#'))].index
        adata = adata[ind_cells]
    return adata

## ===  Greenleaf AML study PBMC GSE139369
def add_AMLPBMC_metadata(data_dir, adata):
    '''Loading cell metadata for AML study PBMC
    '''
    metadata = pd.read_csv(data_dir+os.sep+'GSE139369_PBMC_scATACseq/PBMC_healthy_cells.tsv',
            header=0, sep='\t')
    metadata["ID"] = metadata["Group"]+'#'+metadata["Barcode"]
    metadata["cell.type"] = metadata["BioClassification"].str.split(pat='_', expand=True)[1]  ## extract out cell types
    metadata = metadata[metadata["cell.type"] != "Unk"]  ## remove unknown

    common_barcodes = set(metadata["ID"]).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(metadata, left_on="barcode", right_on="ID", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    return adata

def load_AMLPBMC_genescore(data_dir, ind=None, model_name="GeneModel-GB-Exponential-Extend-18"):
    '''Load AML PBMC gene score matrix derived from ArchR
    '''
    matrice_dir = data_dir+os.sep+'GSE139369_PBMC_scATACseq/ArchR_GeneModels'
    ## load adata
    adata = anndata.read_mtx(matrice_dir+os.sep+model_name+'.mtx.gz').T
    ## load cells and genes
    genes = pd.read_csv(matrice_dir+os.sep+model_name+'_genes.tsv', header=None, sep='\t')
    adata.var["genes"] = genes[0].values
    adata.var_names = adata.var["genes"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None
    cells = pd.read_csv(matrice_dir+os.sep+model_name+'_cells.tsv', header=None, sep='\t')
    adata.obs["barcode"] = [x[18:28]+'#'+x.split('#')[1] for x in cells[0].values.tolist()]
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    adata = add_AMLPBMC_metadata(data_dir, adata)
    ## filter by individual
    if ind != "all" and ind is not None:
        ind_cells = adata.obs[adata.obs["Group"].isin(ind.split('#'))].index
        adata = adata[ind_cells]
 
    return adata

def load_AMLPBMC_ATACbins(data_dir, ind=None):
    '''Load immune PBMC ATAC-seq 500bp bins matrix derived from ArchR

    @ind: replicates separated by # separator (same as how the ArchR did for barcode), can be one of PBMC_D10T1, PBMC_D11T1, PBMC_D12T1, PBMC_D12T2, PBMC_D12T3
    '''
    matrice_dir = data_dir+os.sep+'GSE139369_PBMC_scATACseq/matrices'

    reps = ["PBMC_D10T1", "PBMC_D11T1", "PBMC_D12T1", "PBMC_D12T2", "PBMC_D12T3"]  ## 3 individual with 5 replicates
    if ind == "all" or ind is None:
        inds = reps
    else:
        inds = ind.split('#')

    adata_list = []
    for rep in inds:
        rep_adata = anndata.read_mtx(matrice_dir+os.sep+rep+'_tile.mtx').T
        rep_cells = pd.read_csv(matrice_dir+os.sep+rep+'_tile_barcodes.tsv', 
                header=None, sep='\t')
        rep_adata.obs["barcode"] = rep_cells[0].values
        rep_adata.obs_names = rep_adata.obs['barcode']
        rep_adata.obs_names_make_unique(join="-")
        rep_adata.obs.index.name = None
        adata_list.append(rep_adata)

    from anndata import AnnData
    adata = AnnData.concatenate(*adata_list, index_unique=None, batch_categories=inds)
    ## add cell metadata
    adata = add_AMLPBMC_metadata(data_dir, adata)
    return adata

# === load Greenleaf data all together
def load_Greenleaf_genescore(data_dir, ind=None, model_name="GeneModel-GB-Exponential-Extend-18"):
    immunePBMC_inds, AMLPBMC_inds = None, None
    if '+' in ind:
        immunePBMC_inds, AMLPBMC_inds = ind.split('+')
        immunePBMC_inds = immunePBMC_inds.replace('immunePBMC#', '')
        AMLPBMC_inds = AMLPBMC_inds.replace('AMLPBMC#', '')
    else:
        immunePBMC_inds = ind
        immunePBMC_inds = immunePBMC_inds.replace('immunePBMC#', '')

    adata_list = []
    if immunePBMC_inds is not None:
        immunePBMC_adata = load_immunePBMC_genescore(data_dir, ind=immunePBMC_inds,
            model_name=model_name)
        immunePBMC_adata = curate_immunePBMC_celltypes(immunePBMC_adata, celltype_gran=0)
        adata_list.append(immunePBMC_adata)
    if AMLPBMC_inds is not None:
        AMLPBMC_adata = load_AMLPBMC_genescore(data_dir, ind=AMLPBMC_inds,
            model_name=model_name)
        AMLPBMC_adata = curate_AMLPBMC_celltypes(AMLPBMC_adata, celltype_gran=0)
        adata_list.append(AMLPBMC_adata)
    ## extract out common columns
    common_columns = ['barcode', 'Group', 'ID', 'cell.type']
    for adata in adata_list:
        adata_obs = adata.obs[common_columns]
        adata.obs = adata_obs
    adata = anndata.AnnData.concatenate(*adata_list, join="inner")
    del adata_list  ## release space

    return adata

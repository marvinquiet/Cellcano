import os
import anndata
import numpy as np
import pandas as pd

def load_sciCAR_ge(data_dir, gene_file=None):
    '''Loading sciCAR gene expression

    @ gene_file: whether to load pre-calculated correlated genes
    '''
    ## load data as anndata
    adata = anndata.read_mtx(data_dir+os.sep+"sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_gene_count.txt").T

    ## load genes and cells
    genes = pd.read_csv(data_dir+os.sep+"sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_gene.txt",
            header=0, sep=',')
    genes.drop_duplicates(subset = ['gene_short_name'], keep = 'first', inplace = True)  ## remove duplicates
    protein_genes = genes[genes['gene_type'] == 'protein_coding']
    adata = adata[:, protein_genes.index]
    adata.var["gene_symbols"] = protein_genes["gene_short_name"].values
    adata.var_names = adata.var["gene_symbols"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    ## if an input gene file is given
    if gene_file:
        selected_genes = []
        with open(gene_file) as f:
            selected_genes = f.readlines()
        if len(selected_genes) > 0:  ## whether there are items in selected_genes
            selected_genes = [x.strip() for x in selected_genes]
        common_genes = set(adata.var_names).intersection(set(selected_genes))
        adata = adata[:, list(common_genes)]

    cells = pd.read_csv(data_dir+os.sep+"sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt",
            header=0, sep=',')
    drop_NA_cells = cells[cells["cell_name"].notna()]
    adata = adata[drop_NA_cells.index, :]
    adata.obs['barcode'] = drop_NA_cells["sample"].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    ATAC_cells = pd.read_csv(data_dir+os.sep+"sciCAR_GSE117089/GSM3271045_ATAC_mouse_kidney_cell.txt",
            header=0, sep=',')
    common_barcodes = set(adata.obs["barcode"]).intersection(set(ATAC_cells["sample"]))
    adata = adata[list(common_barcodes)]  ## select common barcodes betweeen RNA and ATAC

    adata.obs = adata.obs.merge(cells, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata

def load_sciCAR_SC2P(data_dir, option="Ftest"):
    '''Load sciCAR SC2P binarized dataset
    '''
    if option == "Ftest":
        adata = anndata.read_csv(data_dir+os.sep+'sciCAR_GSE117089/SC2P_Ftest_exprs.csv').T
    if option == "twoGroupT":
        adata = anndata.read_csv(data_dir+os.sep+'sciCAR_GSE117089/SC2P_twoGroupT_exprs.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load cell metadata
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata


def load_sciCAR_GA(data_dir):
    '''Loading sciCAR gene activity data
    '''
    ## load data as anndata
    adata = anndata.read_mtx(data_dir+os.sep+'sciCAR_GSE117089/sciCAR_cicero_gene_activities.mtx').T

    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/sciCAR_cicero_gene_activities_genes.tsv',
            header=None, sep='\t')
    adata.var['gene_symbols'] = genes[0].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    cells = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/sciCAR_cicero_gene_activities_barcodes.tsv',
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")

    ## load cell information
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None
    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)

    return adata


def load_sciCAR_promoter(data_dir, gene_file=None, binarize=False):
    '''Loading sciCAR promoter counts

    @ gene_file: whether to select out certain genes
    @ binarize: whether to binarize the promoter counts
    '''
    ## load data as anndata
    if binarize:
        adata = anndata.read_mtx(data_dir+os.sep+'sciCAR_GSE117089/promoter_binarized_counts.mtx').T
    else:
        adata = anndata.read_mtx(data_dir+os.sep+'sciCAR_GSE117089/promoter_counts.mtx').T
    ## load genes
    genes = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/promoter_genes.tsv',
            header=None, sep='\t')
    adata.var['gene_symbols'] = genes[0].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    ## filter genes if given
    if gene_file:
        selected_genes = []
        with open(gene_file) as f:
            selected_genes = f.readlines()
        if len(selected_genes) > 0:  ## whether there are items in selected_genes
            selected_genes = [x.strip() for x in selected_genes]
        common_genes = set(adata.var_names).intersection(set(selected_genes))
        adata = adata[:, list(common_genes)]

    ## load cells
    cells = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/promoter_barcodes.tsv',
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")

    ## load cell information
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None
    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)

    return adata

def load_sciCAR_peaks(data_dir):
    '''Loading sciCAR peaks
    '''
    ## load data as anndata
    adata = anndata.read_mtx(data_dir+os.sep+"sciCAR_GSE117089/GSM3271045_ATAC_mouse_kidney_peak_count.txt").T

    ## load cells and peaks
    peaks = pd.read_csv(data_dir+os.sep+"sciCAR_GSE117089/GSM3271045_ATAC_mouse_kidney_peak.txt",
            header=0, sep=',')
    adata.var["peaks"] = peaks["peak"].values
    adata.var_names = adata.var["peaks"]
    adata.var_names_make_unique(join="-")
    adata.var.index.name = None

    cells = pd.read_csv(data_dir+os.sep+"sciCAR_GSE117089/GSM3271045_ATAC_mouse_kidney_cell.txt",
            header=0, sep=',')
    adata.obs['barcode'] = cells["sample"].values
    adata.obs_names = adata.obs['barcode']
    adata.obs_names_make_unique(join="-")
    adata.obs.index.name = None

    ## load cells from RNA-seq to get cell types and get common cells
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata

def load_sciCAR_cusanovich(data_dir):
    '''Loading sciCAR top 20,000 features extracted by Cusanovich2018 by TF-IDF and SVD
    '''
    adata = anndata.read_csv(data_dir+os.sep+'sciCAR_GSE117089/cusanovich_LSI_matrix.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load cell metadata
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata

def load_sciCAR_cisTopic(data_dir):
    '''Loading sciCAR cisTopic results -> cell*7 topics (dealt on Emory RSPH cluster)
    '''
    adata = anndata.read_csv(data_dir+os.sep+'sciCAR_GSE117089/cisTopic_result.csv').T
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load cell metadata
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata

def load_sciCAR_SnapATAC(data_dir, option="raw"):
    '''Loading sciCAR SnapATAC bin-by-cell matrix/ PCs-by-cell matrix

    @option: 
        raw - raw bin-by-cell matrix
        filter - filtered bin-by-cell matrix
        PC - PC-by-cell matrix based on filtered matrix
    '''
    if option == "raw":
        adata = anndata.read_mtx(data_dir+os.sep+'sciCAR_GSE117089/SnapATAC_bin.mtx')
        cells = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/SnapATAC_barcodes.tsv',
                header=None, sep='\t')
        adata.obs['barcode'] = cells[0].values
        adata.obs_names = cells[0]
    elif option == "filter":
        adata = anndata.read_mtx(data_dir+os.sep+'sciCAR_GSE117089/SnapATAC_bin_filtered.mtx')
        cells = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/SnapATAC_barcodes_filtered.tsv',
                header=None, sep='\t')
        adata.obs['barcode'] = cells[0].values
        adata.obs_names = cells[0]
    elif option == "PC":
        adata = anndata.read_csv(data_dir+os.sep+'sciCAR_GSE117089/SnapATAC_50PCs_mat.csv')
        adata.obs['barcode'] = adata.obs.index

    adata.obs_names_make_unique(join="-")
    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")
    ## load cell metadata
    cell_metadata = pd.read_csv(data_dir+os.sep+'sciCAR_GSE117089/GSM3271044_RNA_mouse_kidney_cell.txt',
            header=0, sep=',')
    common_barcodes = set(cell_metadata['sample']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(cell_metadata, left_on="barcode", right_on="sample", how='left')
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell_name': 'cell.type'}, inplace=True)
    return adata



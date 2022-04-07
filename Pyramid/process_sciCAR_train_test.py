import anndata
import random

from ATACseq.preprocess import load_sciCAR_data

def process_sciCAR_ge(data_dir, train_pct=0.8, sample_seed=0, gene_file=None):
    ''' Process sciCAR gene expression

    @ data_dir: where sciCAR gene expression is stored
    @ result_dir: where to store PCA/tSNE/UMAP result
    @ train_pct: training dataset percentage
    '''
    adata = load_sciCAR_data.load_sciCAR_ge(data_dir, gene_file=gene_file)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_sciCAR_SC2P(data_dir, train_pct=0.8, sample_seed=0, option="Ftest"):
    '''Process sciCAR binarized gene expresison data by SC2P estimation

    @option: 
        Ftest: 1000 features selected by F-test
        twoGroupT: 174 features selected by two-group t-test
    '''
    adata = load_sciCAR_data.load_sciCAR_SC2P(data_dir, option=option)

    ## eliminate non-expressed cells/genes
    import scanpy.api as sc
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)

    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata


def process_sciCAR_SC2P_to_promoter(data_dir, option="Ftest"):
    '''Process sciCAR SC2P gene expression to predict binarized promoter counts

    @option:
        Ftest: 1000 features selected by F-test
        twoGroupT: 174 features selected by two-group t-test
    '''
    train_adata = load_sciCAR_data.load_sciCAR_SC2P(data_dir, option=option)
    ## eliminate non-expressed cells/genes
    import os
    import scanpy.api as sc
    sc.pp.filter_cells(train_adata, min_genes=1)
    sc.pp.filter_genes(train_adata, min_cells=1)

    if option == "Ftest":
        gene_file = data_dir+os.sep+'sciCAR_GSE117089/F-test_markergenes.tsv'
    if option == "twoGroupT":
        gene_file=data_dir+os.sep+'sciCAR_GSE117089/twoGroupT_markergenes.tsv'
    test_adata = load_sciCAR_data.load_sciCAR_promoter(data_dir, gene_file=gene_file, binarize=True)
    sc.pp.filter_cells(test_adata, min_genes=1)
    sc.pp.filter_genes(test_adata, min_cells=1)

    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]
    return train_adata, test_adata

def process_sciCAR_GA(data_dir, train_pct=0.8, sample_seed=0):
    ''' Process sciCAR gene activity

    @ data_dir: where sciCAR gene activity is stored
    @ result_dir: where to store PCA/tSNE/UMAP result
    @ train_pct: training dataset percentage
    '''
    adata = load_sciCAR_data.load_sciCAR_GA(data_dir)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_sciCAR_promoter(data_dir, train_pct=0.8, sample_seed=0, 
        gene_file=None, binarize=False):
    '''Process sciCAR promoter counts

    @ gene_file: load data with selected genes
    @ binarize: whether to binarize the promoter counts matrix
    '''
    adata = load_sciCAR_data.load_sciCAR_promoter(data_dir, gene_file=gene_file,
            binarize=binarize)

    ## eliminate non-expressed cells/genes
    import scanpy.api as sc
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)

    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_sciCAR_peaks(data_dir, train_pct=0.8, sample_seed=0):
    '''Process sciCAR peaks
    '''
    adata = load_sciCAR_data.load_sciCAR_peaks(data_dir)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_sciCAR_cusanovich(data_dir, train_pct=0.8, sample_seed=0):
    '''Process sciCAR matrix from cusanovich2018 transformation
    '''
    adata = load_sciCAR_data.load_sciCAR_cusanovich(data_dir)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata


def process_sciCAR_cisTopic(data_dir, train_pct=0.8, sample_seed=0):
    '''Process sciCAR matrix from cisTopic transformation
    '''
    adata = load_sciCAR_data.load_sciCAR_cisTopic(data_dir)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_sciCAR_SnapATAC(data_dir, option="raw", train_pct=0.8, sample_seed=0):
    '''Process sciCAR matrix from SnapATAC bin-by-cell matrix
    '''
    adata = load_sciCAR_data.load_sciCAR_SnapATAC(data_dir, option=option)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

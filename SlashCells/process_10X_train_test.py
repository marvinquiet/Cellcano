import os
import anndata
import random

from ATACseq.preprocess import load_10X_data


def random_cv(adata, train_pct, sample_seed=0):
    '''Random split cross-validation train/test
    '''
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata


def subset_genes(adata, data_dir, gene_opt):
    '''Load specific genes for gene expresion or promoter counts
    @ gene_opt:
        - Pearson: F-test selected marker genes intersect with scATAC-seq
            Pearson correlation genes
        - twoGroupT: F-test selected marker genes intersect with scATAC-seq 
            two group t-test genes
        - Ftest: F-test selected marker genes
    '''
    selected_genes = []
    if gene_opt == "Ftest":
        gene_file = data_dir+os.sep+'10X_Multiome_PBMC10K/F-test_markergenes.tsv'
    if gene_opt == "Pearson":
        gene_file = data_dir+os.sep+'10X_Multiome_PBMC10K/PearsonCor_markergenes.tsv'
    if gene_opt == "twoGroupT":
        gene_file = data_dir+os.sep+'10X_Multiome_PBMC10K/twoGroupT_markergenes.tsv'

    if gene_file:
        with open(gene_file) as f:
            selected_genes = f.readlines()
        if len(selected_genes) > 0:  ## whether there are items in selected_genes
            selected_genes = [x.strip() for x in selected_genes]
            common_genes = set(adata.var_names).intersection(set(selected_genes))
            adata = adata[:, list(common_genes)]

    return adata

def process_10X_ArchRgenescore(data_dir, train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC matrix from ArchR gene scores
    '''
    #adata = load_10X_data.load_10XPBMC_ArchRgenescore(data_dir)
    adata = load_10X_data.load_10XPBMC_genescore(data_dir)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_SnapATAC(data_dir, option="raw", train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC matrix from SnapATAC bin-by-cell matrix, or PC

    @option:
        - raw: bin-by-cell matrix
        - PC: SVD PCs
    '''
    adata = load_10X_data.load_10XPBMC_SnapATAC(data_dir, option=option)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_ATACbin(data_dir, train_pct=0.8, sample_seed=0):
    '''Process 10X ATACbin from ArchR bin-by-cell matrix

    '''
    adata = load_10X_data.load_10XPBMC_ATACbins(data_dir)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_peak(data_dir, train_pct=0.8, sample_seed=0):
    '''Process 10XPBMC matrix from peak-by-cell matrix
    '''
    adata = load_10X_data.load_10XPBMC_peak(data_dir)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_promoter(data_dir, gene_opt=None, binarize=True, train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC matrix from promoter-by-cell matrix
    @binarize: whether to binarize the dataset
    '''
    adata = load_10X_data.load_10XPBMC_promoter(data_dir, binarize=binarize)
    if gene_opt:
        adata = subset_genes(adata, data_dir, gene_opt)

    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_ge(data_dir, gene_opt=None, train_pct=0.8, sample_seed=0):
    '''Process 10XPBMC matrix from gene expresson
    '''
    adata = load_10X_data.load_10XPBMC_ge(data_dir)
    if gene_opt:
        adata = subset_genes(adata, data_dir, gene_opt)

    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_ge_SC2P(data_dir, gene_opt=None, train_pct=0.8, sample_seed=0):
    '''Process 10XPBMC matrix from gene expresson
    '''
    adata = load_10X_data.load_10XPBMC_ge_SC2P(data_dir)
    if gene_opt:
        adata = subset_genes(adata, data_dir, gene_opt)

    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_SC2P_to_binpromoter(data_dir, gene_opt=None):
    '''Use trained models on SC2P to predict binarized promoter
    '''
    train_adata = load_10X_data.load_10XPBMC_ge_SC2P(data_dir)
    test_adata = load_10X_data.load_10XPBMC_promoter(data_dir, binarize=True)
    if gene_opt:
        train_adata = subset_genes(train_adata, data_dir, gene_opt)
        test_adata = subset_genes(test_adata, data_dir, gene_opt)
    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]
    return train_adata, test_adata

def process_10X_exprs_to_promoter(data_dir, gene_opt=None):
    '''Use trained models on gene expression to predict promoter counts
    '''
    train_adata = load_10X_data.load_10XPBMC_ge(data_dir)
    test_adata = load_10X_data.load_10XPBMC_promoter(data_dir, binarize=False)
    if gene_opt:
        train_adata = subset_genes(train_adata, data_dir, gene_opt)
        test_adata = subset_genes(test_adata, data_dir, gene_opt)
    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]
    return train_adata, test_adata


def process_10X_ga(data_dir, gene_opt=None, train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC gene activities matrix
    '''
    adata = load_10X_data.load_10XPBMC_ga(data_dir)
    if gene_opt:
        adata = subset_genes(adata, data_dir, gene_opt)

    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_cusanovich(data_dir, train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC cusanovich transformation
    '''
    adata = load_10X_data.load_10XPBMC_Cusanovich(data_dir)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_cisTopic(data_dir, train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC cisTopic transformation
    '''
    adata = load_10X_data.load_10XPBMC_cisTopic(data_dir)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

def process_10X_Seurat_imputed(data_dir, opt="ga", train_pct=0.8, sample_seed=0):
    ### Process 10X PBMC Seurat imputed scATAC-seq
    train_adata = load_10X_data.load_10XPBMC_Seurat_imputed(data_dir)
    if opt == "ga":
        test_adata = load_10X_data.load_10XPBMC_ga(data_dir)
    if opt == "promoter":
        test_adata = load_10X_data.load_10XPBMC_promoter(data_dir, binarize=False)

    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]
 
    return train_adata, test_adata

def process_10X_Seurat_NNweights(data_dir, opt="WKNN", train_pct=0.8, sample_seed=0):
    '''Process 10X PBMC Seurat weights
    '''
    adata = load_10X_data.load_10XPBMC_Seurat_NNweights(data_dir, opt=opt)
    return random_cv(adata, train_pct=train_pct, sample_seed=sample_seed)

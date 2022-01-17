import os
import anndata
import random

from ATACseq.preprocess import load_Greenleaf_data

def random_cv(adata, train_pct, sample_seed=0):
    '''Random split cross-validation train/test
    '''
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

## === Greenleaf immune PBMC
def process_immunePBMC_genescore(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2"):
    '''Process individual data of Greenleaf immune PBMC
    '''
    train_adata = load_Greenleaf_data.load_immunePBMC_genescore(data_dir, ind=ind1)
    test_adata = load_Greenleaf_data.load_immunePBMC_genescore(data_dir, ind=ind2)
    return train_adata, test_adata

def process_immunePBMC_genescore_cv(data_dir, ind="PBMC_Rep2", pct=0.8):
    adata = load_Greenleaf_data.load_immunePBMC_genescore(data_dir, ind=ind)
    return random_cv(adata, train_pct=pct)

def process_immunePBMC_ATACbins(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2"):
    '''Process individual data of Greenleaf immune PBMC
    '''
    train_adata = load_Greenleaf_data.load_immunePBMC_ATACbins(data_dir, ind=ind1)
    test_adata = load_Greenleaf_data.load_immunePBMC_ATACbins(data_dir, ind=ind2)
    return train_adata, test_adata

def process_immunePBMC_ciceroGA(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2"):
    '''Process individual data of Greenleaf immune PBMC
    '''
    train_adata = load_Greenleaf_data.load_immunePBMC_ciceroGA(data_dir, ind=ind1)
    test_adata = load_Greenleaf_data.load_immunePBMC_ciceroGA(data_dir, ind=ind2)
    return train_adata, test_adata

def process_immunePBMC_SnapATAC(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2",
        celltype_gran=1, harmony=False):
    '''Process individual data of immune PBMC SnapATAC
    @harmony: whether to use batch effect removed data or not
    @ind1: always PBMC_Rep1
    @ind2: always PBMC_Rep2
    '''
    train_adata = load_Greenleaf_data.load_immunePBMC_SnapATAC(data_dir, ind=ind1, harmony=harmony)
    test_adata = load_Greenleaf_data.load_immunePBMC_SnapATAC(data_dir, ind=ind2, harmony=harmony)
    return train_adata, test_adata

## === Greeleaf AML study PBMC
def process_AMLPBMC_genescore(data_dir, ind1="PBMC_D10T1", ind2="PBMC_D11T1"):
    '''Process individual data of Greenleaf AML PBMC
    '''
    train_adata = load_Greenleaf_data.load_AMLPBMC_genescore(data_dir, ind=ind1)
    test_adata = load_Greenleaf_data.load_AMLPBMC_genescore(data_dir, ind=ind2)
    return train_adata, test_adata

def process_AMLPBMC_genescore_cv(data_dir, ind="PBMC_D11T1", pct=0.8):
    '''Process individual data of Greenleaf AML PBMC
    '''
    adata = load_Greenleaf_data.load_AMLPBMC_genescore(data_dir, ind=ind)
    return random_cv(adata, train_pct=pct)

def process_AMLPBMC_ATACbins(data_dir, ind1="PBMC_D10T1", ind2="PBMC_D11T1"):
    '''Process individual data of Greenleaf AML PBMC
    '''
    train_adata = load_Greenleaf_data.load_AMLPBMC_ATACbins(data_dir, ind=ind1)
    test_adata = load_Greenleaf_data.load_AMLPBMC_ATACbins(data_dir, ind=ind2)
    return train_adata, test_adata

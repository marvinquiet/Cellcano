import os
import anndata

from typing import TypeVar
A = TypeVar('anndata')  ## generic for anndata

def _is_valid_file(filepath: str) -> bool:
    '''
    whether the file path exists
    ---
    Input
        - filepath: a path to the file
    ---
    Output
        - boolean indicator
    '''
    if not os.path.exists(filepath):
        ## TODO: add error
        pass

    return True

def _data_loader(mtx_file: str, metadata_file: str, celltype_col: str, 
        sep: str = ',', index_col: int = 0) -> A:
    '''
    Load gene score matrix and metadata information
    ---
    Input:
        - mtx_file: gene score matrix in COO format
        - metadata_file: metadata information in either csv, txt, or tsv file
        - celltype_col: which column can be used to extract cell type information
        - sep, index_col: parameters for reading metadata file in pd.read_csv
    ---
    Output:
        - an anndata object
    '''
    pass

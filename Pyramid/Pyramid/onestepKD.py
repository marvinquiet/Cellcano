'''
Method related utility pipelines
'''
import os, sys 
import logging
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy

## import my package
from utils import _utils

## get the logger
logger = logging.getLogger(__name__)

def load_train(inputfile):
    '''Load training data
    ---
    Input:
        - inputfile: input from argument parser
    ---
    Output:
        - adata: gene score anndata object
    '''
    logger.info("Loading data... \n This may take a while depending on your data size..")
    if '.csv' in inputfile:
        adata = _utils._csv_data_loader(inputfile)
    else:
        adata = _utils._COOmtx_data_loader(inputfile)
    return adata


def train_onestepKD(args):
    '''Train one step KD model
        1. Load train data and metadata
        2. Feature selection
        3. Log-Norm and scale data
        4. Train model and save
    ---
    Input:
        - args: user's input arguments
    '''

    '''Two step prediction
    - Select low entropy cells from each cell type
    - Use low entropy cells as reference to predict the rest target cells
    '''
    print('enter train_onestepKD, lolololo')
    logger.info('Enter train_onestepKD')
    teacher_MLP_DIMS = _utils.Teacher_DIMS if args.teacher_ns is None else args.teacher_ns
    student_MLP_DIMS = _utils.Student_DIMS if args.student_ns is None else args.student_ns

    ## load input data
    train_adata = load_train(args.input)
    metadata = _utils._metadata_loader(args.metadata)

    common_cells = set(train_adata.obs_names).intersection(set(metadata.index))
    print("%d common cells found between input data and metadata." % len(common_cells))
    if len(common_cells) < 100:
        print("There are too few cells. Pyramid might not be accurate.")
    train_adata = train_adata[list(common_cells)]
    train_adata.obs = train_adata.obs.merge(metadata, 
            left_on="barcode", right_index=True, how='left')

    ## Feature selection
    if args.fs == "noFS":
        print("Pyramid will not perform feature selection.\n")
        num_features = train_adata.shape[0]
    else:
        num_features = args.num_features
        if num_features < train_adata.shape[0]:
            print("Number of features is larger than data. Pyramid will not perform feature selection.\n")
            num_features = train_adata.shape[0]

    ## preprocess data and select features
    train_adata = _utils._process_adata(train_adata)
    print("Data shape after processing: %d cells X %d genes" % (train_adata.shape[0], train_adata.shape[1]))
    train_adata = _utils._select_feature(train_adata, tmp_dir=args.output_dir,
            fs_method=args.fs, num_features=num_features)
    train_adata = _utils._scale_data(train_adata)
    _utils._visualize_data(train_adata, args.output_dir, prefix=args.prefix)
    _utils._save_adata(train_adata, args.output_dir, prefix=args.prefix)

    ## get x_train, y_train
    x_train = _utils._extract_data(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()

    ## train a KD model
    teacher = _utils._init_MLP(x_train, y_train, dims=_utils.Teacher_DIMS, 
            seed=_utils.RANDOM_SEED)
    teacher.compile()
    student = initialize_MLP(x_train, y_train, dims=_utils.Student_DIMS,
            seed=_utils.RANDOM_SEED)
    distiller = _utils._run_distiller(x_train, y_train, 
            student_model=student.model,
            teacher_model=teacher.model)
    distller.student.save_weights(args.output_dir+args.prefix+'KD.model')


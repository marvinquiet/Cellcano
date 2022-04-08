'''
Functions related to train models
'''
import os, sys 
import logging

import anndata
from sklearn.preprocessing import OneHotEncoder

## import my package
from Pyramid.utils import _utils

## get the logger
logger = logging.getLogger(__name__)

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
    teacher_MLP_DIMS = _utils.Teacher_DIMS if args.teacher_ns is None else args.teacher_ns
    student_MLP_DIMS = _utils.Student_DIMS if args.student_ns is None else args.student_ns


    if args.anndata is not None:
        logger.info("Load pre-processed anndata h5ad file..")
        train_adata = anndata.read_h5ad(args.anndata)
    else:
        if args.input is None or args.metadata is None:
            sys.exit("Please make sure that both gene score matrix and metadata are provided!")

        ## load input data
        logger.info("Loading data... \n This may take a while depending on your data size..")
        if '.csv' in inputfile:
            train_adata = _utils._csv_data_loader(args.input)
        else:
            train_adata = _utils._COOmtx_data_loader(args.input)
        metadata = _utils._metadata_loader(args.metadata)

        common_cells = set(train_adata.obs_names).intersection(set(metadata.index))
        logger.info("%d common cells found between input data and metadata." % len(common_cells))
        if len(common_cells) < 100:
            logger.warning("There are too few cells. Pyramid might not be accurate.")
        train_adata = train_adata[list(common_cells)]
        train_adata.obs = train_adata.obs.merge(metadata, 
                left_on="barcode", right_index=True, how='left')

        ## preprocess data and select features
        train_adata = _utils._process_adata(train_adata)
        logger.info("Data shape after processing: %d cells X %d genes" % (train_adata.shape[0], train_adata.shape[1]))
        train_adata = _utils._select_feature(train_adata, tmp_dir=args.output_dir,
                fs_method=args.fs, num_features=args.num_features)
        train_adata = _utils._scale_data(train_adata)
        _utils._visualize_data(train_adata, args.output_dir, prefix=args.prefix)
        _utils._save_adata(train_adata, args.output_dir, prefix=args.prefix)

    ## get x_train, y_train
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    logger.debug("Categories information: ", enc.categories_)

    ## train a KD model
    teacher = _utils._init_MLP(x_train, y_train, dims=_utils.Teacher_DIMS, 
            seed=_utils.RANDOM_SEED)
    teacher.compile()
    student = _utils._init_MLP(x_train, y_train, dims=_utils.Student_DIMS,
            seed=_utils.RANDOM_SEED)
    distiller = _utils._run_distiller(x_train, y_train, 
            student_model=student.model,
            teacher_model=teacher.model)
    distiller.student.save(args.output_dir+args.prefix+'KD_model')


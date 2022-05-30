'''
Functions related to predict cell types
'''
import os, sys 
import logging

import tensorflow as tf
import numpy as np
import pandas as pd

## import my package
from Cellcano.utils import _utils

## get the logger
logger = logging.getLogger(__name__)

def predict(args):
    model = tf.keras.models.load_model(args.trained_model)
    feature_file = args.trained_model+os.sep+'features.txt'
    encoder_file = args.trained_model+os.sep+'onehot_encoder.txt'
    if not os.path.exists(feature_file) or not os.path.exists(encoder_file):
        sys.exit("Feature file or encoder mapping does not exist! Please check your tained model was trained successfully.")

    features = pd.read_csv(feature_file, sep='\t', header=0, index_col=0)
    encoders = {}
    with open(encoder_file) as f:
        for line in f:
            line_info = line.strip().split(':')
            encoders[int(line_info[0])] = line_info[1]

    ## load input data
    logger.info("Loading data... \n This may take a while depending on your data size..")
    if '.csv' in args.input:
        test_adata = _utils._csv_data_loader(args.input)
    else:
        test_adata = _utils._COOmtx_data_loader(args.input)
    ## process test adata
    test_adata = _utils._process_adata(test_adata, process_type='test')

    ## fill in the data with the same order of features
    feature_idx = []
    find_cnt = 0
    for feature in features.index:
        find_flag = False
        for test_idx, gene in enumerate(test_adata.var_names):
            if gene == feature:
                feature_idx.append(test_idx)
                find_flag = True
                find_cnt += 1
                break
        if not find_flag:
            feature_idx.append(-1)

    if find_cnt < 0.7*len(features):
        logger.warning("The common feature space between reference dataset and target dataset is too few with %d genes.\n This will result in inaccurate prediction." % find_cnt)
    else:
        logger.info("Common feature space between reference and target: %d genes" % find_cnt)

    if -1 not in feature_idx:
        test_adata = test_adata[:, feature_idx]
        logger.info("Data shape after processing: %d cells X %d genes"  % (test_adata.shape[0], test_adata.shape[1]))
    else:
        logger.error("Some features in the reference dataset are not found in the target dataset.")
        sys.exit(1)
        ## TODO: or maybe I should fill in 0 or average of the profiles for those missing genes

    if test_adata.shape[0] >= 1000:
        ## center scale data by test data -> using feature information from test data and do two-step
        test_adata = _utils._scale_data(test_adata)
        test_data_mat = _utils._extract_adata(test_adata)
    else:
        ## scale data by train data mu/std
        test_data_mat = _utils._extract_adata(test_adata)
        test_adata.var['mean'] = np.mean(test_data_mat, axis=0).reshape(-1, 1)
        test_adata.var['std'] = np.std(test_data_mat, axis=0).reshape(-1, 1)
        test_data_mat = (test_data_mat - np.array(features['mean']))/np.array(features['std'])

    y_pred = tf.nn.softmax(model.predict(test_data_mat)).numpy()
    pred_celltypes = _utils._prob_to_label(y_pred, encoders)
    test_adata.obs[_utils.PredCelltype_COLUMN] = pred_celltypes

    if args.tworound:
        ## if less than 1000 in test data
        if test_adata.shape[0] < 1000:
            logger.warning("Your input cell is less than 1000. For performance, we will not perform two-round strategy on your data.")
        ## when cell number is large enough
        else:
            firstround_COLUMN = 'firstround_' + _utils.PredCelltype_COLUMN
            test_adata.obs[firstround_COLUMN] = pred_celltypes
            entropy = [-np.nansum(y_pred[i]*np.log(y_pred[i])) for i in range(y_pred.shape[0])]
            test_adata.obs['entropy'] = entropy
            test_adata = _utils._select_confident_cells(
                    test_adata, celltype_col=firstround_COLUMN)

            low_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'low')].tolist()
            high_entropy_cells = test_adata.obs_names[np.where(test_adata.obs['entropy_status'] == 'high')].tolist()
            test_ref_adata = test_adata[low_entropy_cells]
            test_tgt_adata = test_adata[high_entropy_cells]

            x_tgt_train = _utils._extract_adata(test_ref_adata)
            y_tgt_train = _utils._label_to_onehot(test_ref_adata.obs.loc[low_entropy_cells, firstround_COLUMN].tolist(),
                    encoders=encoders)
            x_tgt_test = _utils._extract_adata(test_tgt_adata)

            ## teahcer/studenmt model on original celltype label
            teacher = _utils._init_MLP(x_tgt_train, y_tgt_train, dims=_utils.Teacher_DIMS,
                    seed=_utils.RANDOM_SEED)
            teacher.compile()
            teacher.fit(x_tgt_train, y_tgt_train, batch_size=_utils.BATCH_SIZE)
            ## student model -> actually same model, just used the concept of distillation
            student = _utils._init_MLP(x_tgt_train, y_tgt_train, dims=_utils.Student_DIMS, 
                    seed=_utils.RANDOM_SEED)
            # Initialize and compile distiller
            distiller = _utils._run_distiller(x_tgt_train, y_tgt_train, 
                    student_model=student.model,
                    teacher_model=teacher.model)
            y_pred_tgt = tf.nn.softmax(distiller.student.predict(x_tgt_test)).numpy()

            pred_celltypes = _utils._prob_to_label(y_pred_tgt, encoders)
            test_adata.obs.loc[high_entropy_cells, _utils.PredCelltype_COLUMN] = pred_celltypes
 
    test_adata.obs.to_csv(args.output_dir+os.sep+args.prefix+'celltypes.csv')


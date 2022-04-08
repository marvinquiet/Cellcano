'''
Functions related to predict cell types
'''
import os, sys 
import logging

## import my package
from Pyramid.utils import _utils

## get the logger
logger = logging.getLogger(__name__)

def predict(args):
    model = tf.kearas.models.load_model(ars.trained_model)

    ## load input data
    logger.info("Loading data... \n This may take a while depending on your data size..")
    if '.csv' in inputfile:
        test_adata = _utils._csv_data_loader(args.input)
    else:
        test_adata = _utils._COOmtx_data_loader(args.input)

    test_adata = _utils._process_adata(test_adata)
    logger.info("Data shape after processing: %d cells X %d genes"  % (test_adata.shape[0], test_adata.shape[1]))
 


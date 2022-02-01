# TensorFlow and tf.keras
import tensorflow.compat.v1 as tf
from tensorflow.python.ops import array_ops, init_ops

from tensorflow import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Concatenate

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

class MLP(object):
    def __init__(self, dims):
        self.dims = dims
        self.model = None
        self.input_shape = None
        self.n_classes = None
        self.random_seed = 0 ## for reproducibility

    def init_MLP(self, dropout_rate=0.1, optimizer='adam'):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN
        model = Sequential()
        for i in range(len(self.dims)):
            if 0 == i:
                model.add(Dense(self.dims[i], input_shape=self.input_shape, 
                    activation=tf.nn.relu,
                    kernel_initializer=dense_kernel_init))
            else:
                model.add(Dense(self.dims[i], activation=tf.nn.relu,
                    kernel_initializer=dense_kernel_init))
            model.add(Dropout(rate=dropout_rate, input_shape=(self.dims[i],), seed=self.random_state))
        model.add(Dense(self.n_classes, activation=tf.nn.softmax,
            kernel_initializer=dense_kernel_init))

        model.compile(loss=tf.keras.losses.CategoricalCrossentropy(), 
                metrics=["accuracy"], ## show training accuracy
                optimizer=optimizer)
        self.model = model
    
    def fit(self, x_train, y_train, batch_size=16, max_epochs=1000, sample_weight=None):
        ## init MLP model
        self.init_MLP()

        ## add callback with 5 steps no improvement
        callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, callbacks=[callback], verbose=2, 
                sample_weight=sample_weight) # without cross validation

    def predict(self, x_test):
        x_pred = self.model.predict(x_test)
        return x_pred




if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"
    from ATACseq.preprocess import load_Greenleaf_data
    from ATACseq.preprocess.process_train_test_data import *
    from ATACseq.preprocess.process_Greenleaf_train_test import *
    from ATACseq import method_utils, dataloading_utils
    import os, time

    ## set arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    args = parser.parse_args()

    args.n_features = 1000
    args.select_on = None
    args.select_method = None

    ## load data
    result_dir = "/home/wma36/gpu/celltyping_methodTrials/results/result_aggregated_model"
    os.makedirs(result_dir, exist_ok=True)

    ATACbin_train_adata, ATACbin_test_adata = \
            process_immunePBMC_ATACbins(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2")
    ATACbin_train_adata = load_Greenleaf_data.curate_immunePBMC_celltypes(ATACbin_train_adata, celltype_gran=0)
    ATACbin_test_adata = load_Greenleaf_data.curate_immunePBMC_celltypes(ATACbin_test_adata, celltype_gran=0)

    GS_train_adata, GS_test_adata = \
            process_immunePBMC_genescore(data_dir, ind1="PBMC_Rep1", ind2="PBMC_Rep2")
    GS_train_adata = load_Greenleaf_data.curate_immunePBMC_celltypes(GS_train_adata, celltype_gran=0)
    GS_test_adata = load_Greenleaf_data.curate_immunePBMC_celltypes(GS_test_adata, celltype_gran=0)

    ## preprocess data
    ATACbin_train_adata, ATACbin_test_adata = \
            dataloading_utils.process_loaded_data(
                    ATACbin_train_adata, ATACbin_test_adata, result_dir, args=args, 
                    save_data=False, ATACseq=True)
    GS_train_adata, GS_test_adata = \
            dataloading_utils.process_loaded_data(
                    GS_train_adata, GS_test_adata, result_dir, args=args,
                    preprocess=True)

    ## find common cells after preprocessing
    common_train_cells = set(ATACbin_train_adata.obs_names).intersection(GS_train_adata.obs_names)
    ATACbin_train_adata = ATACbin_train_adata[list(common_train_cells)]
    GS_train_adata = GS_train_adata[list(common_train_cells)]
    common_test_cells = set(ATACbin_test_adata.obs_names).intersection(GS_test_adata.obs_names)
    ATACbin_test_adata = ATACbin_test_adata[list(common_test_cells)]
    GS_test_adata = GS_test_adata[list(common_test_cells)]

    celltype_cols = "cell.type"

    ## OneHotEncoding the celltypes
    from sklearn.preprocessing import OneHotEncoder
    enc = OneHotEncoder(handle_unknown='ignore')
    ATACbin_train = np.array(ATACbin_train_adata.X)
    GS_train = np.array(GS_train_adata.X)
    y_train = enc.fit_transform(ATACbin_train_adata.obs[[celltype_cols]]).toarray()

    ATACbin_test = np.array(ATACbin_test_adata.X)
    GS_test = np.array(GS_test_adata.X)
    y_test = enc.transform(ATACbin_test_adata.obs[[celltype_cols]]).toarray()

     ## Hyperparameters for network
    dims = [128, 64, 16]
    batch_size = 128
    ### --- run MLP
    print("\n\n=== MLP\n")
    start = time.time()
    cmlp = CombinedMLP(dims)
    cmlp.fit(ATACbin_train, GS_train, y_train=y_train)
    y_pred = cmlp.predict(ATACbin_test, GS_test)
    end = time.time()
    print("\n\n=== Run time:", end-start)

    n_clusters = len(set(GS_train_adata.obs[celltype_cols]))
    method_utils.analyze_prediction(y_pred.argmax(1), y_test.argmax(1), enc, GS_test_adata, 
            n_clusters, end-start, celltype_cols=celltype_cols, 
            result_dir=result_dir, prefix="test_CombinedMLP")


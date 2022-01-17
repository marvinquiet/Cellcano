# TensorFlow and tf.keras
import tensorflow.compat.v1 as tf
from tensorflow.python.ops import array_ops, init_ops

from tensorflow import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Concatenate

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

#print(tf.__version__)

class MLP(object):
    def __init__(self, dims):
        self.dims = dims
        self.model = None
        self.input_shape = None
        self.n_classes = None
        self.random_seed = 0 ## for reproducibility

    def focal_loss(self, y_true, y_pred, gamma=0, alpha=1):
        '''Adapted from tensorflow addons, but also refer to other links
        #https://github.com/fizyr/keras-retinanet/blob/424671f71da40845e987713544493f4e88f4ea68/keras_retinanet/losses.py
        #https://github.com/ailias/Focal-Loss-implement-on-Tensorflow/blob/master/focal_loss.py
        #https://github.com/umbertogriffo/focal-loss-keras/blob/master/src/loss_function/losses.py
        https://github.com/tensorflow/addons/blob/11aad761bfdc5cd4ef29b29e678e0d081f846fa3/tensorflow_addons/losses/focal_loss.py
        https://github.com/artemmavrin/focal-loss/blob/5b2ca68/src/focal_loss/_categorical_focal_loss.py

        .. math::
            L(y, \hat{\mathbf{p}})
                = -\left(1 - \hat{p}_y\right)^\gamma \log(\hat{p}_y)

        Focal loss for imbalanced classification
        
        @y_true: target data
        @y_pred: predicted data
        @gamma: modulating factor, change the focus from well-classified to hard classes
        @alpha: balancing factor, either 1 or None (inverse frequency, or inverse document frequency)
        '''
        if None == alpha:
            alpha = self.weights

        alpha = np.array(alpha, dtype=np.float32)

        epsilon = keras.backend.epsilon()
        y_pred = keras.backend.clip(y_pred, epsilon, 1.-epsilon)
        ce = -y_true * keras.backend.log(y_pred)

        loss = alpha * keras.backend.pow(1-y_pred, gamma) * ce
        return keras.backend.mean(keras.backend.sum(loss, axis=-1))

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
        #model.compile(loss=[self.focal_loss], 
                metrics=["accuracy"], ## show training accuracy
                optimizer=optimizer)
        self.model = model
    
    def fit(self, x_train, y_train, batch_size=16, max_epochs=1000):
        ## calculate weights for each sample using inverse frequency
        self.weights = inverse_variance_weighting(y_train.argmax(1))
        print(self.weights)

        ## init MLP model
        self.init_MLP()

        ## add callback with 5 steps no improvement
        callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, callbacks=[callback], verbose=2) # without cross validation

    def predict(self, x_test):
        x_pred = self.model.predict(x_test)
        return x_pred


class CombinedMLP(object):
    def __init__(self, dims):
        ## initialize two MLP with different dimensions
        self.dims = dims

    def init_CombinedMLP(self):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN
        ATACbin_model = Sequential()
        ATACbin_model.add(Dense(self.dims[0], input_shape=self.input_ATACbin_shape,
            activation=tf.nn.relu, kernel_initializer=dense_kernel_init))
        ATACbin_model.add(Dropout(rate=self.dropout_rate, input_shape=(self.dims[0],), seed=self.random_state))

        GS_model = Sequential()
        GS_model.add(Dense(self.dims[0], input_shape=self.input_GS_shape,
            activation=tf.nn.relu, kernel_initializer=dense_kernel_init))
        GS_model.add(Dropout(rate=self.dropout_rate, input_shape=(self.dims[0],), seed=self.random_state))

        merged_layers = Concatenate()([ATACbin_model.output, GS_model.output])  ## concatenate to a tensor
        for i in range(1, len(self.dims)):
            if 1 == i:
                out = Dense(self.dims[i], activation=tf.nn.relu,
                        kernel_initializer=dense_kernel_init)(merged_layers)
            else:
                out = Dense(self.dims[i], activation=tf.nn.relu,
                         kernel_initializer=dense_kernel_init)(out)
            out = Dropout(rate=self.dropout_rate, input_shape=(self.dims[i],), seed=self.random_state)(out)
        out = Dense(self.n_classes, activation=tf.nn.softmax,
            kernel_initializer=dense_kernel_init)(out)
 
        model = Model([ATACbin_model.input, GS_model.input], out)
        model.compile(loss=tf.keras.losses.CategoricalCrossentropy(),
                metrics=["accuracy"], optimizer=self.optimizer)
        return model

    def fit(self, ATACbin_train, GS_train, y_train, batch_size=16, max_epochs=1000, dropout_rate=0.1, 
            optimizer="adam", random_state=0):
        ## get input shape and initialize parameters
        self.input_ATACbin_shape = (ATACbin_train.shape[1], )
        self.input_GS_shape = (GS_train.shape[1], )
        self.n_classes = len(set(y_train.argmax(1)))

        ## calculate weights for each sample using inverse frequency
        #self.weights = inverse_variance_weighting(y_train.argmax(1))
        #print(self.weights)

        self.dropout_rate = dropout_rate
        self.optimizer = optimizer
        self.random_state = random_state

        ## init MLP model
        self.model = self.init_CombinedMLP()
        print(self.model.summary())

        ## add callback with 5 steps no improvement
        callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        #self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
        #        validation_split=0.2, callbacks=[callback], verbose=2)

        self.model.fit([ATACbin_train, GS_train], y = y_train, 
                epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, callbacks=[callback], verbose=2) # without cross validation

    def predict(self, ATACbin_test, GS_test):
        x_pred = self.model.predict([ATACbin_test, GS_test])
        return x_pred


def inverse_frequency(logits):
    ''' The inverse frequency of labels

    @logits: a list of labels
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)
    freq = counts/sum(counts)  ## get frequency
    inv_freq = 1/freq   ## get inverse frequency

    return inv_freq/sum(inv_freq) ## normalize to 1


def inverse_document_frequency(logits):
    '''Use inverse document frequency with smooth
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)

    freq = np.log(sum(counts)/(1+counts)) + 1
    return freq/sum(freq)  ## normalize to 1

def inverse_variance_weighting(logits):
    '''Inverse variance weighting, using square root
    '''
    unique, counts = np.unique(logits, return_counts=True)
    print(unique, counts)
    freq = counts/sum(counts)
    inv_freq = np.sqrt(1/freq)  ## square root inverse frequency

    return inv_freq/sum(inv_freq) ## normalize to 1

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


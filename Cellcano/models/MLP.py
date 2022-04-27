import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation

class MLP(object):
    def __init__(self, dims):
        self.dims = dims
        self.model = None
        self.input_shape = None
        self.n_classes = None
        self.random_seed = 0 ## for reproducibility

    def init_MLP_model(self, dropout_rate=0.5, 
            optimizer=keras.optimizers.Adam()):
        dense_kernel_init = keras.initializers.TruncatedNormal(mean=0, stddev=0.1, seed=self.random_state) ## same as GEDFN
        model = Sequential()
        model.add(keras.Input(shape=self.input_shape))
        for i in range(len(self.dims)):
            model.add(Dropout(rate=dropout_rate, seed=self.random_state, name="dropout_"+str(i)))
            model.add(Dense(self.dims[i], kernel_initializer=dense_kernel_init, name="dense_"+str(i)))
            model.add(Activation('relu', name="act_"+str(i)))
        model.add(Dense(self.n_classes, kernel_initializer=dense_kernel_init, name="dense_"+str(i+1)))
        self.model = model

    def fit(self, x_train, y_train, batch_size=16, max_epochs=100, 
            sample_weight=None, class_weight=None):
        ## add callback with 5 steps no improvement
        callback = keras.callbacks.EarlyStopping(monitor='loss', patience=5)
        self.model.fit(x_train, y_train, epochs=max_epochs, batch_size=batch_size, 
                validation_split=0.0, callbacks=[callback], verbose=2, 
                sample_weight=sample_weight, class_weight=class_weight) # without cross validation

    def compile(self, optimizer='adam'):
        self.model.compile(
                loss=keras.losses.CategoricalCrossentropy(from_logits=True),
                metrics=["accuracy"], ## show training accuracy,
                optimizer=optimizer)

    def predict(self, x_test):
        ## with softmax
        x_pred = tf.nn.softmax(self.model.predict(x_test)).numpy()
        return x_pred


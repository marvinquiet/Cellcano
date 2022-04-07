## from Keras Knowledge Distillation: https://keras.io/examples/vision/knowledge_distillation/
## for predicting scATAC-seq using another scATAC-seq
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation

class Distiller(keras.Model):
    def __init__(self, student, teacher):
        super(Distiller, self).__init__()
        self.teacher = teacher
        self.student = student

    def compile(
        self,
        optimizer,
        metrics,
        student_loss_fn,
        distillation_loss_fn,
        alpha=0.1,
        temperature=3,
    ):
        """ Configure the distiller.

        Args:
            optimizer: Keras optimizer for the student weights
            metrics: Keras metrics for evaluation
            student_loss_fn: Loss function of difference between student
                predictions and ground-truth
            distillation_loss_fn: Loss function of difference between soft
                student predictions and soft teacher predictions
            alpha: weight to student_loss_fn and 1-alpha to distillation_loss_fn
            temperature: Temperature for softening probability distributions.
                Larger temperature gives softer distributions.
        """
        super(Distiller, self).compile(optimizer=optimizer, metrics=metrics)
        self.student_loss_fn = student_loss_fn
        self.distillation_loss_fn = distillation_loss_fn
        self.alpha = alpha
        self.temperature = temperature

    def train_step(self, data):
        # Unpack data
        x, y = data

        # Forward pass of teacher
        teacher_predictions = self.teacher(x, training=False)

        with tf.GradientTape() as tape:
            # Forward pass of student
            student_predictions = self.student(x, training=True)

            # Compute losses
            student_loss = self.student_loss_fn(y, student_predictions)
            distillation_loss = self.distillation_loss_fn(
                tf.nn.softmax(teacher_predictions / self.temperature, axis=1),
                tf.nn.softmax(student_predictions / self.temperature, axis=1),
            )
            loss = self.alpha * student_loss + (1 - self.alpha) * distillation_loss

        # Compute gradients
        trainable_vars = self.student.trainable_variables
        gradients = tape.gradient(loss, trainable_vars)

        # Update weights
        self.optimizer.apply_gradients(zip(gradients, trainable_vars))

        # Update the metrics configured in `compile()`.
        self.compiled_metrics.update_state(y, student_predictions)

        # Return a dict of performance
        results = {m.name: m.result() for m in self.metrics}
        results.update(
            {"student_loss": student_loss, "distillation_loss": distillation_loss}
        )
        return results

    def test_step(self, data):
        # Unpack the data
        x, y = data

        # Compute predictions
        y_prediction = self.student(x, training=False)

        # Calculate the loss
        student_loss = self.student_loss_fn(y, y_prediction)

        # Update the metrics.
        self.compiled_metrics.update_state(y, y_prediction)

        # Return a dict of performance
        results = {m.name: m.result() for m in self.metrics}
        results.update({"student_loss": student_loss})
        return results


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


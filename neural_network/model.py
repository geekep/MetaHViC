import tensorflow as tf

from keras import layers as kl
from keras import Model, Layer

from MetaHViC.neural_network import igloo, stack_igloo


class OneHotLayer(Layer):
    def call(self, x, word_size):
        return tf.one_hot(x, depth=4**word_size+1, axis=-1)


def create_encoder(input_shape, word_size):
    inputs = kl.Input(shape=(input_shape,), dtype="int64")
    embedded = OneHotLayer()(x=inputs, word_size=word_size)
    outputs = igloo.IGLOO1D_Block(
        embedded,
        nb_patches=2100,
        nb_filters_conv1d=128,
        nb_stacks=3,
        conv1d_kernel=6,
        pooling_size=8,
        dropout_rate=0.2,
        l2_reg=1e-3,
        transformer_style=True,
    )
    outputs = kl.Dense(512)(outputs)
    outputs = kl.BatchNormalization()(outputs)
    outputs = kl.Activation("relu")(outputs)
    return Model(inputs=inputs, outputs=outputs)

def create_classifier(input_shape, word_size, n_classes):
    encoder = create_encoder(input_shape=input_shape, word_size=word_size)
    for layer in encoder.layers:
        layer.trainable = False
    inputs = kl.Input(shape=(input_shape,))
    features = encoder(inputs)
    features = kl.Dense(512)(features)
    features = kl.BatchNormalization()(features)
    features = kl.Activation("relu")(features)
    features = kl.Dropout(0.2)(features)
    outputs = kl.Dense(n_classes, activation="softmax")(features)
    return Model(inputs=inputs, outputs=outputs)


def create_host_encoder(input_shape, word_size):
    inputs = kl.Input(shape=(input_shape,), dtype="int64")
    embedded = OneHotLayer()(x=inputs, word_size=word_size)
    outputs = stack_igloo.IGLOO1D_Block(
        incoming_layer=embedded,
        incoming_proj=0,
        nb_filters_conv1d=128,
        conv1d_kernel=6,
        padding_style="causal",
        spatial_dropout=True,
        nb_stacks=2,
        patch_size=4,
        nb_patches=2100,
        dropout_rate=0.2,
        l2_reg=1e-3,
        transformer_style=True,
        pooling_size=8,
    )
    outputs = kl.Dense(512)(outputs)
    outputs = kl.BatchNormalization()(outputs)
    outputs = kl.Activation("relu")(outputs)
    return Model(inputs=inputs, outputs=outputs)

def create_host_classifier(input_shape, word_size, n_classes):
    encoder = create_host_encoder(input_shape=input_shape, word_size=word_size)
    for layer in encoder.layers:
        layer.trainable = False
    inputs = kl.Input(shape=(input_shape,))
    features = encoder(inputs)
    features = kl.Dense(512)(features)
    features = kl.BatchNormalization()(features)
    features = kl.Activation("relu")(features)
    features = kl.Dropout(0.2)(features)
    outputs = kl.Dense(n_classes, activation="softmax")(features)
    return Model(inputs=inputs, outputs=outputs)

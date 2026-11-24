import os
import re
import gctf
import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow.keras.callbacks import TensorBoard
from keras import Model
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from keras.layers import (
    Activation,
    BatchNormalization,
    Dense,
    Input,
)
from keras.optimizers import Adam


AUTOTUNE = tf.data.experimental.AUTOTUNE
AUTOTUNE = 64
strategy = tf.distribute.MirroredStrategy()
# BATCH_SIZE = 128 * strategy.num_replicas_in_sync
BATCH_SIZE = 64 * strategy.num_replicas_in_sync

FRAGMENT_LENGTH = 6_000
WORD_SIZE = 4
INPUT_SHAPE = FRAGMENT_LENGTH - WORD_SIZE + 1
N_CLASSES = 5
embedding_dim = 128 # set embedding_dim < 512


def parse_tfrecord(example):
    feature_description = {
        "sequence": tf.io.FixedLenFeature(shape=[INPUT_SHAPE], dtype=tf.int64),
        "class": tf.io.FixedLenFeature(shape=[], dtype=tf.int64),
        "sample_weight": tf.io.FixedLenFeature(shape=[], dtype=tf.float32),
    }
    return tf.io.parse_single_example(example, feature_description)

def prepare_sample(features):
    return features["sequence"], features["class"], features["sample_weight"]

def get_dataset(filenames, batch_size):
    return (
        tf.data.TFRecordDataset(filenames, num_parallel_reads=AUTOTUNE)
        .map(map_func=parse_tfrecord, num_parallel_calls=AUTOTUNE)
        .map(map_func=prepare_sample, num_parallel_calls=AUTOTUNE)
        .shuffle(buffer_size=batch_size * 10)
        .batch(batch_size=batch_size, drop_remainder=True)
        .prefetch(buffer_size=AUTOTUNE)
    )

def debug_print(sequence, label, sample_weight):
        tf.print("\nClasses per batch:", label, summarize=-1)
        return sequence, label, sample_weight

def get_dataset_shuffle(filenames, batch_size):
    dataset = tf.data.Dataset.list_files(filenames, shuffle=True)
    dataset = dataset.interleave(
        lambda fname: tf.data.TFRecordDataset(fname),
        cycle_length=AUTOTUNE,
        num_parallel_calls=AUTOTUNE,
        deterministic=False,
    )
    return (
        dataset
        .map(map_func=parse_tfrecord, num_parallel_calls=AUTOTUNE)
        .map(map_func=prepare_sample, num_parallel_calls=AUTOTUNE)
        .shuffle(buffer_size=batch_size * 10)
        .repeat()
        .batch(batch_size=batch_size, drop_remainder=True)
        .map(map_func=debug_print, num_parallel_calls=AUTOTUNE)
        .prefetch(buffer_size=AUTOTUNE)
    )


class SupervisedContrastiveLoss(tf.keras.losses.Loss):
    def __init__(self, temperature=1, name=None):
        super(SupervisedContrastiveLoss, self).__init__(name=name)
        self.temperature = temperature

    def __call__(self, labels, feature_vectors, sample_weight=None):
        feature_vectors_normalized = tf.math.l2_normalize(feature_vectors, axis=1)
        logits = tf.divide(
            tf.matmul(feature_vectors_normalized, tf.transpose(feature_vectors_normalized)),
            self.temperature,
        )
        return tfa.losses.npairs_loss(tf.squeeze(labels), logits)

class SupervisedContrastiveWithCenterLoss(tf.keras.losses.Loss):
    def __init__(self,
                 num_classes,
                 embedding_dim,
                 temperature=1,
                 center_weight=0.5,
                #  center_momentum=0.05,
                 name=None):
        super(SupervisedContrastiveWithCenterLoss, self).__init__(name=name)
        self.temperature = temperature
        self.center_weight = center_weight
        self.centers = tf.Variable(
            initial_value=tf.keras.initializers.GlorotUniform()(shape=(num_classes, embedding_dim)),
            trainable=True,
            name='centers'
        )
        # self.center_momentum = center_momentum
        # self.centers = tf.Variable(
        #     initial_value=tf.keras.initializers.GlorotUniform()(shape=(num_classes, embedding_dim)),
        #     trainable=False,
        #     name="centers"
        # )
    # def _update_centers(self, z, labels):
    #         for k in tf.unique(labels)[0]:
    #             mask = tf.equal(labels, k)
    #             z_k = tf.boolean_mask(z, mask)
    #             if tf.shape(z_k)[0] > 0:
    #                 mean_z = tf.reduce_mean(z_k, axis=0)
    #                 new_center = ((1.0 - self.center_momentum) * self.centers[k] + self.center_momentum * mean_z)
    #                 self.centers[k].assign(new_center)

    def __call__(self, labels, feature_vectors, sample_weight=None):
        feature_vectors_normalized = tf.math.l2_normalize(feature_vectors, axis=1)
        logits = tf.divide(
            tf.matmul(feature_vectors_normalized, tf.transpose(feature_vectors_normalized)),
            self.temperature,
        )
        supcon_loss = tfa.losses.npairs_loss(tf.squeeze(labels), logits)

        batch_labels = tf.cast(tf.squeeze(labels), tf.int32)
        centers_batch = tf.gather(tf.math.l2_normalize(self.centers, axis=1), batch_labels)
        center_loss = tf.reduce_mean(tf.square(feature_vectors_normalized - centers_batch))
        # center_loss = 1.0 - tf.reduce_mean(tf.reduce_sum(feature_vectors_normalized * centers_batch, axis=1))
        # center_loss = tf.reduce_mean(tf.square(feature_vectors_normalized - tf.stop_gradient(centers_batch)))
        # self._update_centers(feature_vectors_normalized, labels)
        
        return supcon_loss + self.center_weight * center_loss


def add_projection_head(encoder):
    inputs = Input(shape=INPUT_SHAPE)
    features = encoder([inputs])
    outputs = Dense(embedding_dim)(features)
    outputs = Activation("relu")(outputs)
    return Model(inputs=[inputs], outputs=outputs)

def create_encoder_with_projection_head(encoder):
    encoder_with_projection_head = add_projection_head(encoder)
    opt = Adam(
        learning_rate=1e-3,
        decay=1e-5,
        clipnorm=1.0,
    )
    opt.get_gradients = gctf.centralized_gradients_for_optimizer(opt)
    encoder_with_projection_head.compile(
        optimizer=opt,
        # loss=SupervisedContrastiveLoss(temperature=0.1),
        loss=SupervisedContrastiveWithCenterLoss(num_classes=N_CLASSES, embedding_dim=embedding_dim, temperature=0.1, center_weight=0.2),
        weighted_metrics=[],
    )
    return encoder_with_projection_head

from stack_igloo import IGLOO1D_Block
def create_encoder():
    inputs = Input(shape=INPUT_SHAPE, dtype="int64")
    onehot = tf.one_hot(inputs, depth=4**WORD_SIZE+1, axis=-1)
    outputs = IGLOO1D_Block(
        incoming_layer=onehot,
        incoming_proj=0,
        nb_filters_conv1d=128,
        conv1d_kernel=6,
        padding_style="causal",
        spatial_dropout=True,
        nb_stacks=3,
        patch_size=4,
        nb_patches=2100,
        dropout_rate=0.2,
        l2_reg=1e-3,
        transformer_style=True,
        pooling_size=8,
    )
    outputs = Dense(512)(outputs)
    outputs = BatchNormalization()(outputs)
    outputs = Activation("relu")(outputs)
    return Model(inputs=inputs, outputs=outputs)

def create_full_model():
    encoder = create_encoder()
    return create_encoder_with_projection_head(encoder)


# checkpoint_dir_encoder = os.path.join("/home/xuzhen/DeepMetaClass/training_checkpoints",
#                                       "encoder_with_projection_head_training_IMG_M_CenterLoss")
checkpoint_dir_encoder = os.path.join("/data/DeepMicroClass/training_checkpoints",
                                      "encoder_with_projection_head_training_IMG_M_CenterLoss")
os.makedirs(checkpoint_dir_encoder, exist_ok=True)
checkpoint_path_encoder = os.path.join(checkpoint_dir_encoder, "encoder-epoch-{epoch:03d}-loss-{loss:.5f}.ckpt")
log_file = os.path.join(checkpoint_dir_encoder, "encoder_training_log_epoch.csv")

cp_callback = ModelCheckpoint(
    filepath=checkpoint_path_encoder,
    save_weights_only=True,
    save_best_only=True,
    save_freq='epoch',
    monitor="val_loss",
    mode="min",
    verbose=1,
)
es_callback = EarlyStopping(
    restore_best_weights=True,
    monitor="val_loss",
    mode='min',
    patience=10,
    min_delta=0,
    verbose=1,
)
csv_logger = CSVLogger(filename=log_file, append=True)
tensorboard = TensorBoard(log_dir=os.path.join(checkpoint_dir_encoder, "logs_batch"),
                          histogram_freq=1,
                          write_graph=True,
                          write_images=False,
                          update_freq=500,
                          profile_batch=(2, 32),
                          embeddings_freq=1,
                          embeddings_metadata=os.path.join(checkpoint_dir_encoder, "logs_batch", "metadata"))

with strategy.scope():
    encoder = create_full_model()
    latest_encoder = tf.train.latest_checkpoint(checkpoint_dir=checkpoint_dir_encoder)
    last_epoch = 0
    if latest_encoder:
        last_epoch = int(re.search(pattern=r"encoder-epoch-(\d+)-loss-([\d.]+)\.ckpt", string=latest_encoder).group(1))
        encoder.load_weights(filepath=latest_encoder)

# tfrecords_dir = "/home/xuzhen/DeepMetaClass/data/tfrecords_length_IMG_M"
tfrecords_dir = "/data/DeepMicroClass/data/tfrecords_length_IMG_M"
train_tfrecords_dir = [f"{tfrecords_dir}/train/*.tfrec",]
train_filenames = sum([tf.io.gfile.glob(dir) for dir in train_tfrecords_dir], [])
val_tfrecords_dir = [f"{tfrecords_dir}/validate/*.tfrec",]
val_filenames = sum([tf.io.gfile.glob(dir) for dir in val_tfrecords_dir], [])

history = encoder.fit(
    x=get_dataset_shuffle(filenames=train_filenames, batch_size=BATCH_SIZE),
    validation_data=get_dataset_shuffle(filenames=val_filenames, batch_size=BATCH_SIZE),
    initial_epoch=last_epoch,
    epochs=150,
    steps_per_epoch=1024,
    validation_steps=128,
    callbacks=[cp_callback, es_callback, csv_logger, tensorboard],
)
encoder.save_weights(os.path.join(checkpoint_dir_encoder, "final_encoder_weights.h5"))

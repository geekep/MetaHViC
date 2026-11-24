import os
import sys
import re
import gctf
import numpy as np
import pandas as pd
from functools import partial
from typing import Optional, Union
from focal_loss import SparseCategoricalFocalLoss
import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow.keras.callbacks import TensorBoard
from keras import Model
from keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from keras.layers import (
    Activation,
    BatchNormalization,
    Dense,
    Dropout,
    Input,
)
from keras.optimizers import Adam

import sequence
sys.path.append(os.path.dirname("DeepMicroClass"))
from DeepMicroClass import constants

AUTOTUNE = tf.data.experimental.AUTOTUNE
strategy = tf.distribute.MirroredStrategy()
BATCH_SIZE = 64 * strategy.num_replicas_in_sync
steps_train_per_epoch = 512
validation_steps = 64

FRAGMENT_LENGTH = 6_000
MIN_LENGTH = 2_500
WORD_SIZE = 4
INPUT_SHAPE = FRAGMENT_LENGTH - WORD_SIZE + 1
N_CLASSES = 5
embedding_dim = 128 # set embedding_dim < 512


def parse_tfrecord(example):
    feature_description = {
        "sequence": tf.io.FixedLenFeature(shape=[INPUT_SHAPE], dtype=tf.int64),
        "class": tf.io.FixedLenFeature(shape=[], dtype=tf.int64),
    }
    return tf.io.parse_single_example(example, feature_description)

def prepare_sample(features):
    return features["sequence"], features["class"]

def get_dataset(filenames, batch_size):
    return (
        tf.data.TFRecordDataset(filenames, num_parallel_reads=AUTOTUNE)
        .map(map_func=parse_tfrecord, num_parallel_calls=AUTOTUNE)
        .map(map_func=prepare_sample, num_parallel_calls=AUTOTUNE)
        .shuffle(buffer_size=batch_size * 10)
        .batch(batch_size=batch_size)
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
        .batch(batch_size=batch_size, drop_remainder=True)
        .repeat()
        .map(map_func=debug_print, num_parallel_calls=AUTOTUNE)
        .prefetch(buffer_size=AUTOTUNE)
    )

def read_fasta_custom_collate(genome_id, label, onehot_dir_list:list, contig_len:Optional[Union[list, int]], word_size:int):
    genome_id = genome_id.numpy().decode('utf-8') if hasattr(genome_id, 'numpy') else str(genome_id)
    for onehot_dir in onehot_dir_list:
        if os.path.exists(os.path.join(onehot_dir, f'{genome_id}.fasta')):
            genome = sequence.read_fasta(filepath=os.path.join(onehot_dir, f'{genome_id}.fasta'), strip_n=True)
            genome = next(genome)
            if contig_len is None:
                contig_len = np.random.choice(a=constants.POSSIBLE_LEN)
            elif isinstance(contig_len, list):
                contig_len = np.random.choice(a=contig_len)
            if contig_len < len(genome):
                start = np.random.randint(low=0, high=len(genome)-contig_len)
                genome = genome[start : start+contig_len]
            contig = sequence.tokenize_dna(seq=genome.seq_ascii.ljust(contig_len, b"N"), word_size=word_size)
            return np.array(contig, dtype=np.int64), np.array([label], dtype=np.int64)
    if contig_len is None:
        contig_len = np.random.choice(a=constants.POSSIBLE_LEN)
    elif isinstance(contig_len, list):
        contig_len = np.random.choice(a=contig_len)
    return np.zeros(contig_len, dtype=np.int64), np.array([-1], dtype=np.int64)

def map_with_shape(genomes_id, labels, onehot_dir_list:list, contig_len:Optional[Union[list, int]], word_size:int):
    contigs, labels = tf.py_function(func=partial(read_fasta_custom_collate, onehot_dir_list=onehot_dir_list, contig_len=contig_len, word_size=word_size),
                                     inp=[genomes_id, labels], Tout=[tf.int64, tf.int64])
    contigs.set_shape([None])
    labels.set_shape([1])
    return contigs, labels

def debug_print(sequence, label):
    tf.print("\nClasses per batch:", label, summarize=-1)
    return sequence, label

def get_dataset_undersampling(sequence_info: str, batch_size: int, onehot_dir_list:list, contig_len:Optional[Union[list, int]], word_size:int):
    sequence_info_df = pd.read_table(sequence_info, sep=',', header=0)
    all_genomes_id = sequence_info_df['genome'].values.astype(str)
    all_class_ = sequence_info_df['class'].values.astype(np.int32)
    unique_classes = np.unique(all_class_)
    class_datasets = []
    for class_ in unique_classes:
        mask = (all_class_ == class_)
        genomes_id = all_genomes_id[mask]
        labels = all_class_[mask]
        class_dataset = tf.data.Dataset.from_tensor_slices(tensors=(genomes_id, labels)).shuffle(buffer_size=len(genomes_id), reshuffle_each_iteration=True).repeat()
        class_datasets.append(class_dataset)
    dataset = tf.data.experimental.sample_from_datasets(datasets=class_datasets, weights=[1/len(unique_classes)]*len(unique_classes))
    return (
        dataset
        .map(partial(map_with_shape, onehot_dir_list=onehot_dir_list, contig_len=contig_len, word_size=word_size), num_parallel_calls=AUTOTUNE)
        .filter(lambda contigs, labels: tf.not_equal(labels[0], -1)) 
        .shuffle(buffer_size=batch_size * 10)
        .batch(batch_size=batch_size, drop_remainder=True)
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
    def __init__(self, num_classes, embedding_dim, temperature=1, center_weight=0.5, name=None):
        super(SupervisedContrastiveWithCenterLoss, self).__init__(name=name)
        self.temperature = temperature
        self.center_weight = center_weight
        self.num_classes = num_classes
        self.centers = self.add_weight(
            name='centers',
            shape=(num_classes, embedding_dim),
            initializer='glorot_uniform',
            trainable=True
        )

    def __call__(self, labels, feature_vectors, sample_weight=None):
        feature_vectors_normalized = tf.math.l2_normalize(feature_vectors, axis=1)
        logits = tf.divide(
            tf.matmul(feature_vectors_normalized, tf.transpose(feature_vectors_normalized)),
            self.temperature,
        )
        supcon_loss = tfa.losses.npairs_loss(tf.squeeze(labels), logits)
        batch_labels = tf.cast(tf.squeeze(labels), tf.int32)
        centers_batch = tf.gather(self.centers, batch_labels)
        center_loss = tf.reduce_mean(tf.square(feature_vectors - centers_batch))
        return supcon_loss + self.center_weight * center_loss

class SparseF1Score(tfa.metrics.F1Score):
    def update_state(self, y_true, y_pred, sample_weight=None):
        y_true = tf.one_hot(tf.cast(y_true, tf.int64), depth=self.num_classes)
        return super().update_state(y_true, y_pred, sample_weight)

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
        loss=SupervisedContrastiveLoss(temperature=0.1),
        # loss=SupervisedContrastiveWithCenterLoss(num_classes=N_CLASSES, embedding_dim=embedding_dim, temperature=0.1, center_weight=0.2),
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
        nb_stacks=1,
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


def create_classifier(encoder, trainable=False):
    for layer in encoder.layers: layer.trainable = trainable
    inputs = Input(shape=INPUT_SHAPE)
    features = encoder(inputs)
    features = Dense(512)(features)
    features = BatchNormalization()(features)
    features = Activation("relu")(features)
    features = Dropout(0.1)(features)
    outputs = Dense(N_CLASSES, activation="softmax")(features)
    classifier = Model(inputs=[inputs], outputs=outputs)
    opt = Adam(learning_rate=7.5e-3, decay=1e-5,)
    opt.get_gradients = gctf.centralized_gradients_for_optimizer(opt)
    classifier.compile(
        optimizer=opt,
        loss=SparseCategoricalFocalLoss(gamma=1),
        weighted_metrics=[tf.keras.metrics.SparseCategoricalAccuracy(name="accuracy"),
                        #   SparseF1Score(num_classes=N_CLASSES, average='macro'),
                        #   tf.keras.metrics.Precision(class_id=1),
                        #   tf.keras.metrics.Recall(class_id=1),
                          ],
    )
    return classifier


# checkpoint_dir_encoder = os.path.join("/home/xuzhen/DeepMetaClass/training_checkpoints",
#                                       "encoder_with_projection_head_training_undersampling")
checkpoint_dir_encoder = os.path.join("/data/DeepMicroClass/training_checkpoints",
                                      "encoder_with_projection_head_training_undersampling_nb_stacks_1_fungi")
os.makedirs(checkpoint_dir_encoder, exist_ok=True)
checkpoint_path_encoder = os.path.join(checkpoint_dir_encoder, "encoder-epoch-{epoch:03d}-loss-{loss:.5f}.ckpt")

# checkpoint_dir_classifier = os.path.join("/home/xuzhen/DeepMetaClass/training_checkpoints",
#                                          "classifier_undersampling")
checkpoint_dir_classifier = os.path.join("/data/DeepMicroClass/training_checkpoints",
                                         "classifier_undersampling_nb_stacks_1_fungi")
os.makedirs(checkpoint_dir_classifier, exist_ok=True)
checkpoint_path_classifier = os.path.join(checkpoint_dir_classifier, "classifier-epoch-{epoch:03d}-loss-{loss:.5f}.ckpt")
log_file = os.path.join(checkpoint_dir_classifier, "classifier_training_log_epoch.csv")

cp_callback = ModelCheckpoint(
    filepath=checkpoint_path_classifier,
    save_weights_only=True,
    save_best_only=True,
    save_freq='epoch',
    monitor="val_loss",
    mode="min",
    verbose=1,
)
es_callback = EarlyStopping(
    restore_best_weights=True,
    monitor="val_accuracy",
    mode='max',
    patience=10,
    min_delta=0,
    verbose=1,
)
csv_logger = CSVLogger(filename=log_file, append=True)
tensorboard = TensorBoard(log_dir=os.path.join(checkpoint_dir_classifier, "logs_batch"),
                          histogram_freq=1,
                          write_graph=True,
                          write_images=False,
                          update_freq=500,
                          profile_batch=(2, 32),
                          embeddings_freq=1,
                          embeddings_metadata=os.path.join(checkpoint_dir_classifier, "logs_batch", "metadata"))

with strategy.scope():
    encoder = create_encoder()
    encoder_with_projection_head = create_encoder_with_projection_head(encoder)
    # latest_encoder = tf.train.latest_checkpoint(checkpoint_dir=checkpoint_dir_encoder)
    # encoder_with_projection_head.load_weights(filepath=latest_encoder)
    encoder_with_projection_head.load_weights(filepath=os.path.join(checkpoint_dir_encoder, "final_encoder_weights.h5"))
    classifier = create_classifier(encoder)
    latest_classifier = tf.train.latest_checkpoint(checkpoint_dir=checkpoint_dir_classifier)
    last_epoch = 0
    if latest_classifier:
        last_epoch = int(re.search(pattern=r"classifier-epoch-(\d+)-loss-([\d.]+)\.ckpt", string=latest_classifier).group(1))
        classifier.load_weights(filepath=latest_classifier)

INPUT_PATH = '/data/DeepMicroClass/data'
ONEHOT_PATH_LIST = [os.path.join(INPUT_PATH, f) for f in ['pre20_prokaryote.fa.gz_single_fasta',
                                                          'prok_vir_Train_Validate.fasta.gz_single_fasta',
                                                          'pre20_plasmid.fa.gz_single_fasta',
                                                          'pre20_fungi.fa.gz_single_fasta',
                                                          'euk_vir_Train_Validate.fasta.gz_single_fasta']]
TRAINING_DATASET_CSV = os.path.join(INPUT_PATH, 'train_class_5_fungi.csv')
VAL_DATASET_CSV = os.path.join(INPUT_PATH, 'validate_class_5_fungi.csv')

# train_class_5_weight = pd.read_table(TRAINING_DATASET_CSV, sep=',', header=0)
# classes = np.unique(train_class_5_weight['class'])
# weights = train_class_5_weight['windows_per_genome'].sum() / len(classes) / train_class_5_weight.groupby('class')['windows_per_genome'].sum()
# classes2weights = dict(zip(classes, weights))

history = classifier.fit(
    x=get_dataset_undersampling(sequence_info=TRAINING_DATASET_CSV, batch_size=BATCH_SIZE, onehot_dir_list=ONEHOT_PATH_LIST, contig_len=FRAGMENT_LENGTH, word_size=WORD_SIZE),
    validation_data=get_dataset_undersampling(sequence_info=VAL_DATASET_CSV, batch_size=BATCH_SIZE, onehot_dir_list=ONEHOT_PATH_LIST, contig_len=FRAGMENT_LENGTH, word_size=WORD_SIZE),
    initial_epoch=last_epoch,
    epochs=120,
    steps_per_epoch=steps_train_per_epoch,
    validation_steps=validation_steps,
    callbacks=[cp_callback, es_callback, csv_logger, tensorboard],
)
classifier.save_weights(os.path.join(checkpoint_dir_classifier, "final_classifier_weights.h5"))

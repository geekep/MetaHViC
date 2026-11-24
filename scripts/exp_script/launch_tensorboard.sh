#!/bin/bash

# eval "$(conda shell.bash hook)"
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate gctf

# train and validate set from tfrecords_length_6000
tensorboard --logdir_spec igloo-encoder-6000:/data/DeepMetaClass/training_checkpoints/encoder_with_projection_head_training/logs_batch,igloo-classifier-6000:/data/DeepMetaClass/training_checkpoints/classifier/logs_batch --port 6006
# train and validate set from tfrecords_length_IMG_M
tensorboard --logdir_spec igloo-encoder-IMG-M:/data/DeepMetaClass/training_checkpoints/encoder_with_projection_head_training_IMG_M/logs_batch,igloo-classifier-IMG-M:/data/DeepMetaClass/training_checkpoints/classifier_IMG_M/logs_batch --port 6007
# train and validate set from undersampling
tensorboard --logdir_spec igloo-encoder-6000-undersampling:/data/DeepMetaClass/training_checkpoints/encoder_with_projection_head_training_undersampling/logs_batch,igloo-classifier-6000-undersampling:/data/DeepMetaClass/training_checkpoints/classifier_undersampling/logs_batch --port 6008
# train and validate set from undersampling (including Fungi)
tensorboard --logdir_spec igloo-encoder-6000-undersampling_fungi:/data/DeepMetaClass/training_checkpoints/encoder_with_projection_head_training_undersampling_nb_stacks_1_fungi/logs_batch,igloo-classifier-6000-undersampling_fungi:/data/DeepMetaClass/training_checkpoints/classifier_undersampling_nb_stacks_1_fungi/logs_batch --port 6009

# mu03
tensorboard --logdir_spec igloo-encoder-IMG-M:/home/xuzhen/DeepMetaClass/training_checkpoints/encoder_with_projection_head_training_IMG_M/logs_batch,igloo-classifier-IMG-M:/home/xuzhen/DeepMetaClass/training_checkpoints/classifier_IMG_M/logs_batch --port 6007
ssh -L 6007:localhost:6007 xuzhen@58.206.100.240 -p 22
import os
import sys
import gzip
import argparse
import math
import numpy as np
import pandas as pd
from Bio import SeqIO, Data
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument('--train_input_file', dest='train_input_file', required=True)
parser.add_argument('--validate_input_file', dest='validate_input_file', required=True)
parser.add_argument('--output_dir', dest='output_dir')
parser.add_argument('--class_id', dest='class_id', type=int, choices=[0, 1, 2, 3, 4], required=True)
parser.add_argument('--to_csv_header', dest='to_csv_header', action='store_true')
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.train_input_file))
os.makedirs(args.output_dir, exist_ok=True)

df_train = []
for recoed in SeqIO.parse(gzip.open(args.train_input_file, 'rt'), 'fasta'):
    df_train.append([recoed.id] + [str(args.class_id)])
print(len(df_train))
with open(os.path.join(args.output_dir, 'train_class_5_fungi.csv'), 'at') as out:
    df = pd.DataFrame(data=df_train, columns=['sequence', 'class'])
    if args.to_csv_header:
        df.to_csv(out, sep=',', index=False, header=True)
    else:
        df.to_csv(out, sep=',', index=False, header=False)

df_validate = []
for recoed in SeqIO.parse(gzip.open(args.validate_input_file, 'rt'), 'fasta'):
    df_validate.append([recoed.id] + [str(args.class_id)])
print(len(df_validate))
with open(os.path.join(args.output_dir, 'validate_class_5_fungi.csv'), 'at') as out:
    df = pd.DataFrame(data=df_validate, columns=['sequence', 'class'])
    if args.to_csv_header:
        df.to_csv(out, sep=',', index=False, header=True)
    else:
        df.to_csv(out, sep=',', index=False, header=False)

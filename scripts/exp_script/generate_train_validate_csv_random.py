import os, sys
import gzip
import argparse
import math
import numpy as np
import pandas as pd
from Bio import SeqIO, Data
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument('--input_file', dest='input_file', required=True)
parser.add_argument('--output_dir', dest='output_dir')
parser.add_argument("--validate_ratio", dest="validate_ratio", default=0.2, type=float)
parser.add_argument('--class_id', dest='class_id', type=int, choices=[0, 1, 2, 3, 4], required=True)
parser.add_argument('--to_csv_header', dest='to_csv_header', action='store_true')
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.train_input_file))
os.makedirs(args.output_dir, exist_ok=True)

org_records = SeqIO.parse(gzip.open(args.input_file, 'rt'), "fasta")
org_id = [record.id for record in org_records]
shuffled_index = np.arange(len(org_id))
np.random.shuffle(shuffled_index)
boundary = math.floor(len(org_id) * (1 - float(args.validate_ratio)))

train_id = [org_id[i] for i in shuffled_index[:boundary]]
df_train = [[id, str(args.class_id)] for id in set(train_id)]
print(len(df_train))
with open(os.path.join(args.output_dir, 'train_class_5.csv'), 'a') as out:
    df = pd.DataFrame(data=df_train, columns=['sequence', 'class'])
    if args.to_csv_header:
        df.to_csv(out, sep=',', index=False, header=True)
    else:
        df.to_csv(out, sep=',', index=False, header=False)

val_id = [org_id[i] for i in shuffled_index[boundary:]]
df_validate = [[id, str(args.class_id)] for id in set(val_id)]
print(len(df_validate))
with open(os.path.join(args.output_dir, 'validate_class_5.csv'), 'a') as out:
    df = pd.DataFrame(data=df_validate, columns=['sequence', 'class'])
    if args.to_csv_header:
        df.to_csv(out, sep=',', index=False, header=True)
    else:
        df.to_csv(out, sep=',', index=False, header=False)

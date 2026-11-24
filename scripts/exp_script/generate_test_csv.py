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
parser.add_argument('--test_input_file', dest='test_input_file', required=True)
parser.add_argument('--output_dir', dest='output_dir')
parser.add_argument('--class_id', dest='class_id', type=int, choices=[0, 1, 2, 3, 4], required=True)
parser.add_argument('--to_csv_header', dest='to_csv_header', action='store_true')
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.train_input_file))
os.makedirs(args.output_dir, exist_ok=True)

df_test = []
for recoed in SeqIO.parse(gzip.open(args.test_input_file, 'rt'), 'fasta'):
    df_test.append([recoed.id] + [str(args.class_id)])
print(len(df_test))
with open(os.path.join(args.output_dir, 'test_class_5.csv'), 'at') as out:
    df = pd.DataFrame(data=df_test, columns=['genome', 'class'])
    if args.to_csv_header:
        df.to_csv(out, sep=',', index=False, header=True)
    else:
        df.to_csv(out, sep=',', index=False, header=False)

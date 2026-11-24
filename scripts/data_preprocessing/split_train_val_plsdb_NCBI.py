import argparse
import math
import os
import numpy as np
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--input", dest="input", required=True)
parser.add_argument("--output_dir", dest="output_dir")
parser.add_argument("--train_ratio", dest="train_ratio", default=0.8, type=float)
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.input))

with gzip.open(args.input, 'rt') as handle:
    org_id = [record.id for record in SeqIO.parse(handle, "fasta")]
shuffled_index = np.arange(len(org_id))
np.random.shuffle(shuffled_index)
boundary = math.floor(len(org_id) * args.train_ratio)
train_id = set(org_id[i] for i in shuffled_index[:boundary])
# val_id = set(org_id[i] for i in shuffled_index[boundary:])

with gzip.open(args.input, 'rt') as handle, \
     gzip.open(f"{args.output_dir}/{str(os.path.basename(args.input)).rsplit('.')[0]}_Train.gz", 'at') as train_out, \
     gzip.open(f"{args.output_dir}/{str(os.path.basename(args.input)).rsplit('.')[0]}_Validate.gz", 'at') as val_out:
    
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in train_id:
            SeqIO.write(record, train_out, "fasta")
        else:
            SeqIO.write(record, val_out, "fasta")

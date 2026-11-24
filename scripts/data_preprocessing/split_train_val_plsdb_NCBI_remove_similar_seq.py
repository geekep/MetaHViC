import argparse
import math
import pandas as pd
import os
import numpy as np
import gzip
from Bio import SeqIO
# from ete3 import NCBITaxa


parser = argparse.ArgumentParser()
parser.add_argument("--input", dest="input", required=True)
parser.add_argument("--output_dir", dest="output_dir")
parser.add_argument("--train_ratio", dest="train_ratio", default=0.8, type=float)
parser.add_argument("--dist", dest="dist", default="", type=str)
parser.add_argument("--threshold", dest="threshold", default=0.1, type=float)
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.input))

with gzip.open(args.input, 'rt') as handle:
    org_records = list(SeqIO.parse(handle=handle, format="fasta"))
org_id = [record.id for record in org_records]
shuffled_index = np.arange(len(org_id))
np.random.shuffle(shuffled_index)
boundary = math.floor(len(org_id) * args.train_ratio)
train_id = [org_id[i] for i in shuffled_index[:boundary]]
val_id = [org_id[i] for i in shuffled_index[boundary:]]

dist = pd.read_table(args.dist, index_col=0, delimiter="\t", engine="c", na_filter=False, low_memory=False)
train_val_dist = dist.loc[train_id, val_id]
train_val_overlap = set(train_val_dist.columns[(train_val_dist < args.threshold).sum(axis=0) == 0])
train_id.extend(train_val_overlap)
train_id = set(train_id)
val_id = set(id for id in val_id if id not in train_val_overlap)

train_records = (record for record in org_records if record.id in train_id)
with gzip.open(f"{args.output_dir}/{str(os.path.basename(args.input)).split('.', 1)[0]}_Train.fasta.gz", 'wt') as out:
    SeqIO.write(train_records, out, "fasta")

val_records = (record for record in org_records if record.id in val_id)
with gzip.open(f"{args.output_dir}/{str(os.path.basename(args.input)).split('.', 1)[0]}_Validate.fasta.gz", 'wt') as out:
    SeqIO.write(val_records, out, "fasta")

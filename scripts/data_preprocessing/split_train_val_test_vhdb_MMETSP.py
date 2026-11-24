import argparse
import gzip
import math
from Bio import SeqIO
import pandas as pd
import os
import numpy as np
from random import shuffle
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", dest="input_dir", required=True)
parser.add_argument("--gzip", dest="use_gzip", action='store_true')
parser.add_argument("--output_dir", dest="output_dir")
parser.add_argument("--train_ratio", dest="train_ratio", default=0.6, type=float)
parser.add_argument("--validate_ratio", dest="validate_ratio", default=0.2, type=float)
args = parser.parse_args()

input_list = [fn for fn in os.listdir(args.input_dir) if fn.endswith((".gz", ".fasta", ".fna", ".fa")) and not fn.startswith('.')]
org_records = []
for f in tqdm(input_list):
    if args.use_gzip: input_handle = gzip.open(os.path.join(args.input_dir, f), 'rt')
    else: input_handle = open(os.path.join(args.input_dir, f), "rt")
    org_records.extend(list(SeqIO.parse(input_handle, format="fasta")))
    
org_id = [record.id for record in org_records]
shuffle(org_id)
border1 = math.floor(len(org_id) * args.train_ratio)
border2 = border1 + math.floor(len(org_id) * args.validate_ratio)
train_ids = set(org_id[:border1])
val_ids = set(org_id[border1:border2])
test_ids = set(org_id[border2:])

if not args.output_dir: args.output_dir = os.path.abspath(args.input_dir)
with gzip.open(f"{args.output_dir}/{os.path.basename(args.input_dir)}_Train.fasta.gz", "at") as train_handle, \
     gzip.open(f"{args.output_dir}/{os.path.basename(args.input_dir)}_Validate.fasta.gz", "at") as val_handle, \
     gzip.open(f"{args.output_dir}/{os.path.basename(args.input_dir)}_Test.fasta.gz", "at") as test_handle:
    
    for record in org_records:
        if record.id in train_ids:
            SeqIO.write(record, train_handle, "fasta")
        elif record.id in val_ids:
            SeqIO.write(record, val_handle, "fasta")
        else:
            SeqIO.write(record, test_handle, "fasta")

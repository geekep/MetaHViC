from Bio import SeqIO
import os, sys
import argparse
import gzip
import numpy as np

sys.path.append(os.path.dirname("DeepMicroClass"))
from DeepMicroClass import utils

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', dest='input', required=True)
parser.add_argument('--output_dir', '-o', dest='output_dir')
args = parser.parse_args()

if not args.output_dir: args.output_dir = os.path.dirname(os.path.abspath(args.input))
os.makedirs(args.output_dir, exist_ok=True)
split_dir_fasta = os.path.join(args.output_dir, f'{os.path.basename(os.path.abspath(args.input))}_single_fasta')
os.makedirs(split_dir_fasta, exist_ok=True)
# split_dir_numpy = os.path.join(args.output_dir, f'{os.path.basename(os.path.abspath(args.input))}_single_numpy')
# os.makedirs(split_dir_numpy, exist_ok=True)

for record in SeqIO.parse(gzip.open(args.input, 'rt'), 'fasta'):
    refseqid = record.id.replace('/', '')
    # refseqid = record.description.split(' ')[0].replace('/', '')

    if not os.path.exists(os.path.join(split_dir_fasta, f'{refseqid}.fasta')):
        SeqIO.write(record, os.path.join(split_dir_fasta, f'{refseqid}.fasta'), 'fasta')

    # if not os.path.exists(os.path.join(split_dir_numpy, f'{refseqid}.npy')):
    #     seq = ''.join(str(record.seq))
    #     onehot = utils.seq2onehot(seq, num_classes=5)[:, :4]
    #     np.save(os.path.join(split_dir_numpy, f'{refseqid}'), onehot)

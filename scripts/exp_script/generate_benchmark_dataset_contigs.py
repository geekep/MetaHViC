import os
import gzip
import argparse
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", dest="input_dir", required=True)
parser.add_argument("--output_dir", "-o", dest="output_dir", required=True)
parser.add_argument("--length_max", "-l", dest="length_max", default=6_000, type=int, help="The maximum length of the sampled contigs")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

def sample_seq(seq, length):
    if len(seq) <= length:
        return seq
    else:
        start = np.random.randint(0, len(seq)-length)
        return seq[start : start+length]

for f in tqdm([fn for fn in os.listdir(args.input_dir) if fn.endswith((".gz", ".fasta", ".fna", ".fa")) and not fn.startswith('.')]):
    if f.endswith('.gz'):
        input_generator = SeqIO.parse(gzip.open(os.path.join(args.input_dir, f), 'rt'), 'fasta')
    else:
        input_generator = SeqIO.parse(open(os.path.join(args.input_dir, f), 'rt'), 'fasta')
    output_record = [sample_seq(record, int(args.length_max)) for record in input_generator]
    with gzip.open(os.path.join(args.output_dir, f), 'wt') as out:
        SeqIO.write(output_record, out, 'fasta')

from Bio import SeqIO
from tqdm import tqdm
import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', '-i', dest='input_file', required=True, nargs='+',)
# parser.add_argument('--input_file', '-i', dest='input_file', required=True, action='append',)
# parser.add_argument('--input_file', '-i', dest='input_file', required=True, type=lambda s: s.split(','),)
parser.add_argument("--gzip", "-z", dest="use_gzip", action='store_true')
parser.add_argument('--output_dir', '-o', dest='output_dir', required=True)
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)
input_list = [fn for fn in args.input_file if fn.endswith((".gz", ".fasta", ".fna", ".fa")) and not fn.startswith('.')]

recoed_id = set()
recoed_length = len(recoed_id)
recoed = list()
for _, fasta_file in enumerate(input_list):
    if args.use_gzip: input_handle = gzip.open(fasta_file, "rt")
    else: input_handle = open(fasta_file, "r")
    for seq in tqdm(SeqIO.parse(input_handle, format="fasta"), leave=False):
        recoed_id.add(seq.id.split('.')[0])
        if len(recoed_id) > recoed_length:
            recoed.append(seq)
            recoed_length += 1
    input_handle.close()

with gzip.open(os.path.join(f"{args.output_dir}, merged.fasta.gz")) as out:
    SeqIO.write(recoed, out, "fasta")

import argparse
import gzip
import os
from Bio import SeqIO
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", dest="input_dir")
parser.add_argument("--gzip", "-z", dest="use_gzip", action='store_true')
parser.add_argument("--output_dir", "-o", dest="output_dir")
parser.add_argument("--plsdb_path", "-p", dest="plsdb_path", default='/data/DeepMicroClass/data/plsdb/plsdb.fna.gz', type=str)
args = parser.parse_args()

# Check if input is same as output
assert args.input_dir != args.output_dir, "Input and output directory cannot be the same"
os.makedirs(args.output_dir, exist_ok=True)

with gzip.open(args.plsdb_path) as handle: plasmid_ids = set(seq.id for seq in SeqIO.parse(handle, "fasta"))
input_list = [fn for fn in os.listdir(args.input_dir) if fn.endswith((".gz")) and not fn.startswith('.')]
for fasta_file in tqdm(input_list):
    input_fasta_path = os.path.join(args.input_dir, fasta_file)
    if args.use_gzip: input_handle = gzip.open(input_fasta_path, "rt")
    else: input_handle = open(input_fasta_path, "rt")
    seq_iter = SeqIO.parse(input_handle, format="fasta")
    seq_list = []
    for seq in tqdm(seq_iter, leave=False):
        if "plasmid" not in seq.description.lower() and (seq.id not in plasmid_ids):
            seq_list.append(seq)
    input_handle.close()
    
    output_fasta_path = os.path.join(args.output_dir, fasta_file)
    if args.use_gzip: output_handle = gzip.open(output_fasta_path, "wt")
    else: output_handle = open(output_fasta_path, "wt")
    SeqIO.write(seq_list, output_handle, "fasta")
    output_handle.close()

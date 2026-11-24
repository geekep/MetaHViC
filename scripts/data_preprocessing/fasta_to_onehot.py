import sys
sys.path.append('..')
import numpy as np
import util
import os
from Bio import SeqIO, Data
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', dest='input_dir', required=True)
parser.add_argument('--output_dir', dest='output_dir')
parser.add_argument('--amino_acids', dest='amino_acids', action='store_true')
args = parser.parse_args()

if not args.output_dir: args.output_dir = f'{args.input_dir}_single_onehot'
os.makedirs(args.output_dir, exist_ok=True)

possible_aa = set(Data.IUPACData.protein_letters)
possible_aa.add('*')
possible_aa.add('X')
possible_aa.add('B')
possible_aa.add('Z')
possible_aa.add('J')
possible_aa.add('O')
possible_aa.add('U')
aa_index = {aa: i for i, aa in enumerate(possible_aa)}

input_fasta = os.listdir(args.input_dir)
input_fasta = list(filter(lambda x: x.endswith(('.fna', '.fasta', '.fa', '.gz')), input_fasta))
for fasta in tqdm(input_fasta):
    seq_records = list(SeqIO.parse(os.path.join(args.input_dir, fasta), 'fasta'))

    output_onehot_path = os.path.join(args.output_dir, fasta + '.npy')
    if not os.path.exists(output_onehot_path):
        seq = [str(record.seq) for record in seq_records]
        seq = ''.join(seq)
        onehot = util.seq2onehot(seq, num_classes=5)[:, :4]
        np.save(output_onehot_path, onehot.astype(np.int8))
    
    if args.amino_acids:
        os.makedirs(args.output_dir + '_aa', exist_ok=True)
        output_onehot_aa_path = os.path.join(args.output_dir + '_aa', fasta + '.npy')
        if not os.path.exists(output_onehot_aa_path):
            aa_sequences = [record.translate() for record in seq_records]
            aa_sequences = [str(record.seq) for record in aa_sequences]
            aa_sequences1 = [record[1:].translate() for record in seq_records]
            aa_sequences1 = [str(record.seq) for record in aa_sequences1]
            aa_sequences2 = [record[2:].translate() for record in seq_records]
            aa_sequences2 = [str(record.seq) for record in aa_sequences2]
            aa_sequences = ''.join(aa_sequences + aa_sequences1 + aa_sequences2)
            aa_indices = [aa_index[aa] for aa in aa_sequences]
            aa_onehot = np.zeros((len(aa_sequences), len(aa_indices)), dtype=np.int8)
            aa_onehot[np.arange(len(aa_indices)), aa_indices] = 1
            np.save(output_onehot_aa_path, aa_onehot.astype(np.int8))
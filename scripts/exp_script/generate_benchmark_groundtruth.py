import os, re, gzip
import numpy as np
import pandas as pd
from Bio import SeqIO

FASTA_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_fungi'
results_fn = list(filter(lambda x: x.endswith(('.gz', '.fasta', '.fna', '.fa')) and not x.startswith('.'), os.listdir(FASTA_DIR)))
results_fn.sort()

groundtruth_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_fungi_groundtruth'
os.makedirs(groundtruth_DIR, exist_ok=True)

def construct_target(nums):
    order = [[3], [4], [2], [0], [1]]
    return np.array(sum([order[i] * nums[i] for i in range(len(order))], []))

for f in results_fn:
    nums = [int(n) for n in re.findall(r'\d+', f)[1:]]
    target = construct_target(nums)
    if f.endswith('.gz'):
        ids = [seq.id for seq in SeqIO.parse(gzip.open(os.path.join(FASTA_DIR, f), 'rt'), 'fasta')]
    else:
        ids = [seq.id for seq in SeqIO.parse(open(os.path.join(FASTA_DIR, f), 'rt'), 'fasta')]
    df = pd.DataFrame(data=target, index=ids)
    df.to_pickle(path=os.path.join(groundtruth_DIR, f'{f}.pkl'))

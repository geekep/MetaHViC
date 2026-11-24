import os
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from random import sample
from matplotlib import pyplot as plt

PROK_PATH = "/data/DeepMicroClass/data/clean_sequence_filtered_fasta_0.1/prokaryote.fa.gz"
PROKVIR_PATH = "/data/DeepMicroClass/data/clean_sequence_filtered_fasta_0.1/prokvir.fa.gz"
PLASMID_PATH = "/data/DeepMicroClass/data/clean_sequence_filtered_fasta_0.1/plasmid.fa.gz"
EUK_PATH = "/data/DeepMicroClass/data/clean_sequence_filtered_fasta_0.1/fungi.fa.gz"
EUKVIR_PATH = "/data/DeepMicroClass/data/clean_sequence_filtered_fasta_0.1/eukvir.fa.gz"
PATHS = [PROK_PATH, PROKVIR_PATH, PLASMID_PATH, EUK_PATH, EUKVIR_PATH]
records = [list(record for record in SeqIO.parse(gzip.open(path, 'rt'), 'fasta') if len(record) > 0) for path in PATHS]
print([len(r) for r in records])

TARGET_PATH = "/data/DeepMicroClass/data/sequence_subsampled_1000_class_5_fungi"
os.makedirs(TARGET_PATH, exist_ok=True)
TOTAL_SEQ_NUM = 1000

PROK_EUK_RATIO = np.array([[5, 1],
                           [4, 2],
                           [3, 3],
                           [2, 4],
                           [1, 5]], dtype=float)
PROK_EUK_RATIO /= PROK_EUK_RATIO.sum(axis=1, keepdims=True)

PROK_PROKVIR_PLASMID_RATIO = np.array([
    [4, 1, 1],
    [3, 1, 1],
    [2, 1, 1],
    [1, 1, 1],
], dtype=float)
PROK_PROKVIR_PLASMID_RATIO /= PROK_PROKVIR_PLASMID_RATIO.sum(axis=1, keepdims=True)

EUK_EUKVIR_RATIO = np.array([
    [4, 2],
    [3, 3],
    [2, 4],
    [1, 5],
], dtype=float)
EUK_EUKVIR_RATIO /= EUK_EUKVIR_RATIO.sum(axis=1, keepdims=True)

ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * PROK_PROKVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0)
seq_num = np.round(ratio * TOTAL_SEQ_NUM).astype(int)
print(seq_num)
for i in range(seq_num.shape[0]):
    sampled_records = [sample(population=records[j], k=seq_num[i, j]) for j in range(seq_num.shape[1])]
    output_path = os.path.join(TARGET_PATH,
                               f'{i:02d}_PROK{len(sampled_records[0])}_PROKVIR{len(sampled_records[1])}_PLSMD{len(sampled_records[2])}_EUK{len(sampled_records[3])}_EUKVIR{len(sampled_records[4])}.fa.gz')
    SeqIO.write(sum(sampled_records, []),
                gzip.open(output_path, 'wt'),
                'fasta')


# plot ds_composition
xticks = [f'DS_{i}' for i in range(1, 21)]
ratio_df = pd.DataFrame(ratio, columns=['Prokaryote', 'ProkaryoteVirus', 'Plasmid', 'Fungi', 'EukaryoteVirus'], index=xticks)
ax = ratio_df.plot(kind='bar',
                   stacked=True,
                   width=0.8,
                   color=["#b48a1f", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
                   xlabel="Dataset",
                   ylabel="Percentage (%)")
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
ax.legend(loc='lower left', prop={'size': 10})
plt.xticks(rotation=45)
plt.savefig('plot_figures/plot_benchmark/analysis_results/benchmark_dataset_composition_fungi_bar.pdf', bbox_inches='tight')

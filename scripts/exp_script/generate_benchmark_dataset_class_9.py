import os
import math
import numpy as np
import pandas as pd
from Bio import SeqIO
from random import sample
from matplotlib import pyplot as plt

BACT_PATH = "data/clean_sequence_filtered_fasta_0.1/bacteria.gz"
BACTVIR_PATH = "data/clean_sequence_filtered_fasta_0.1/Bacteria_vir.fa.gz"
ARCH_PATH = "data/clean_sequence_filtered_fasta_0.1/archaea.fa.gz"
ARCHVIR_PATH = "data/clean_sequence_filtered_fasta_0.1/ArchaeaVirus.fa.gz"
PLASMID_PATH = "data/clean_sequence_filtered_fasta_0.1/plasmid.fa.gz"
FUNGI_PATH = "data/clean_sequence_filtered_fasta_0.1/fungi.fa.gz"
FUNGIVIR_PATH = "data/clean_sequence_filtered_fasta_0.1/FungiVirus.fa.gz"
EUK_PATH = "data/clean_sequence_filtered_fasta_0.1/eukaryote.fa.gz"
EUKVIR_PATH = "data/clean_sequence_filtered_fasta_0.1/eukvir.fa.gz"
PATHS = [BACT_PATH, BACTVIR_PATH, ARCH_PATH, ARCHVIR_PATH, PLASMID_PATH, FUNGI_PATH, FUNGIVIR_PATH, EUK_PATH, EUKVIR_PATH]

TARGET_PATH = "data/sequence_subsampled_1000_class_9"
os.makedirs(TARGET_PATH, exist_ok=True)

TOTAL_SEQ_NUM = 1000

PROK_EUK_RATIO = np.array([[9, 1],
                           [7, 3],
                           [5, 5],
                           [3, 7],
                           [1, 9]], dtype=float)
PROK_EUK_RATIO /= PROK_EUK_RATIO.sum(axis=1, keepdims=True)

BACT_BACTVIR_ARCH_ARCHVIR_PLASMID_RATIO = np.array([[9, 1, 9, 0.3, 9],
                                                    [7, 1, 7, 0.3, 7],
                                                    [5, 1, 5, 0.3, 5],
                                                    [3, 1, 3, 0.3, 3]], dtype=float)
BACT_BACTVIR_ARCH_ARCHVIR_PLASMID_RATIO /= BACT_BACTVIR_ARCH_ARCHVIR_PLASMID_RATIO.sum(axis=1, keepdims=True)

FUNGI_FUNGIVIR_EUK_EUKVIR_RATIO = np.array([[9, 1, 9, 1],
                                            [7, 1, 7, 1],
                                            [5, 1, 5, 1],
                                            [3, 1, 3, 1]], dtype=float)
FUNGI_FUNGIVIR_EUK_EUKVIR_RATIO /= FUNGI_FUNGIVIR_EUK_EUKVIR_RATIO.sum(axis=1, keepdims=True)

ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * BACT_BACTVIR_ARCH_ARCHVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * FUNGI_FUNGIVIR_EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0)
seq_num = np.round(ratio * TOTAL_SEQ_NUM).astype(int)
records = [[record for record in SeqIO.parse(path, 'fasta') if len(record) > 5000] for path in PATHS]
for i in range(seq_num.shape[0]):
    sampled_records = [sample(records[j], seq_num[i, j]) for j in range(seq_num.shape[1])]
    output_path = os.path.join(TARGET_PATH, f'{i:02d}_Bacteria-{len(sampled_records[0])}_BacteriaVirus-{len(sampled_records[1])}_Archaea-{len(sampled_records[2])}_ArchaeaVirus-{len(sampled_records[3])}_Plasmid-{len(sampled_records[4])}_Fungi-{len(sampled_records[5])}_FungiVirus-{len(sampled_records[6])}_Eukaryote-{len(sampled_records[7])}_EukaryoteVirus-{len(sampled_records[8])}.fa')
    SeqIO.write(sum(sampled_records, []), output_path, 'fasta')


xticks = [f'DS_{i}' for i in range(1, 21)]
ratio_df = pd.DataFrame(ratio, columns=['Bacteria', 'BacteriaVirus', 'Archaea', 'ArchaeaVirus', 'Plasmid', 'Fungi', 'FungiVirus', 'Eukaryote', 'EukaryoteVirus'], index=xticks)
ax = ratio_df.plot(kind='bar',
                   stacked=True,
                   width=0.8,
                   color=["#b48a1f", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', "#1f35b4", "#b41fa8", "#b44c1f", "#82b41f"],
                   xlabel="Dataset",
                   ylabel="Percentage (%)")
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
ax.legend(loc='lower left', prop={'size': 10})
plt.xticks(rotation=45)
plt.savefig('figures/ds_composition.pdf', bbox_inches='tight')

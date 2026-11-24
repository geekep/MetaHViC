from Bio import SeqIO
import pandas as pd
import os
import gzip
import numpy as np

DB_METADATA_PATH = "data/plsdb/meta/nuccore.csv"
FASTA_PATH = "data/plsdb/plsdb.fna"
PRE20_PATH = "data/plasmid/pre20_plasmid.fa.gz"
os.makedirs(os.path.dirname(os.path.abspath(PRE20_PATH)), exist_ok=True)
POST20_PATH = "data/plasmid/post20_plasmid.fa.gz"
os.makedirs(os.path.dirname(os.path.abspath(POST20_PATH)), exist_ok=True)

metadata = pd.read_csv(DB_METADATA_PATH, sep=',')
dates = pd.to_datetime(metadata['NUCCORE_CreateDate'])
sep_date = np.datetime64('2020-01-01')
pre20_set = set(metadata.loc[dates < sep_date]['NUCCORE_ACC'])
post20_set = set(metadata.loc[dates >= sep_date]['NUCCORE_ACC'])

pre20_records = []
post20_records = []
for record in SeqIO.parse(FASTA_PATH, 'fasta'):
    if record.id in pre20_set: pre20_records.append(record)
    elif record.id in post20_set: post20_records.append(record)

with gzip.open(PRE20_PATH) as handle:
    SeqIO.write(pre20_records, handle, 'fasta')
with gzip.open(POST20_PATH) as handle:
    SeqIO.write(post20_records, handle, 'fasta')
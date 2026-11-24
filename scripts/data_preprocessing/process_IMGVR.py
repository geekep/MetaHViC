from Bio import SeqIO
import pandas as pd
import os
import gzip
# from ete3 import NCBITaxa

IMGVR_FASTA_PATH = "/data/DeepMicroClass/data/IMGVR/IMGVR_all_nucleotides.fna.gz"

DB_METADATA_PATH = "/data/DeepMicroClass/data/IMGVR/IMGVR_all_Host_information.tsv"
metadata = pd.read_table(DB_METADATA_PATH)
metadata.dropna(subset=["Host taxonomy prediction"], inplace=True)
metadata.reset_index(drop=True, inplace=True)

prok_vir = set()
Archaea_vir = set()
Bacteria_vir = set()
euk_vir = set()
for i, row in metadata.iterrows():
    host = row["Host taxonomy prediction"].split(";")[0]
    UVIG_id = row["UVIG"]
    if "Archaea" in host or "Bacteria" in host:
        prok_vir.add(UVIG_id)
        if   "Archaea" in host:  Archaea_vir.add(UVIG_id)
        elif "Bacteria" in host: Bacteria_vir.add(UVIG_id)
    elif "Eukaryota" in host:
        euk_vir.add(UVIG_id)

OUTPUT_Host_PATH = os.path.join(os.path.dirname(os.path.abspath(IMGVR_FASTA_PATH)), "Host")
os.makedirs(OUTPUT_Host_PATH, exist_ok=True)
with gzip.open(IMGVR_FASTA_PATH, 'rt') as input_handle, \
     gzip.open(os.path.join(OUTPUT_Host_PATH, "prok_vir.fasta.gz"), "at") as prok_vir_handle, \
     gzip.open(os.path.join(OUTPUT_Host_PATH, "Archaea_vir.fasta.gz"), "at") as Archaea_vir_handle, \
     gzip.open(os.path.join(OUTPUT_Host_PATH, "Bacteria_vir.fasta.gz"), "at") as Bacteria_vir_handle, \
     gzip.open(os.path.join(OUTPUT_Host_PATH, "euk_vir.fasta.gz"), "at") as euk_vir_handle:

    for record in SeqIO.parse(input_handle, "fasta"):
        UVIG = record.id.split("|")[0]
        if UVIG in prok_vir:
            SeqIO.write(record, prok_vir_handle, "fasta")
            if UVIG in Archaea_vir:
                SeqIO.write(record, Archaea_vir_handle, "fasta")
            elif UVIG in Bacteria_vir:
                SeqIO.write(record, Bacteria_vir_handle, "fasta")
        elif UVIG in euk_vir:
            SeqIO.write(record, euk_vir_handle, "fasta")


# DB_METADATA_PATH = "/data/DeepMicroClass/data/IMGVR/IMGVR_all_Sequence_information.tsv"
# metadata = pd.read_table(DB_METADATA_PATH)
# metadata.dropna(subset=["Ecosystem classification"], inplace=True)
# metadata.reset_index(drop=True, inplace=True)

# conds = ['Animal', 'Aquatic_Sediment', 'Plants']
# conds_sets = [
#     set(metadata[metadata["Ecosystem classification"].str.contains(conds[0])]["UVIG"]),
#     set(metadata[metadata["Ecosystem classification"].str.contains(conds[1])]["UVIG"]),
#     set(metadata[metadata["Ecosystem classification"].str.contains(conds[2])]["UVIG"]),
# ]

# OUTPUT_Ecosystem_PATH = os.path.join(os.path.dirname(os.path.abspath(IMGVR_FASTA_PATH)), "Ecosystem")
# os.makedirs(OUTPUT_Ecosystem_PATH, exist_ok=True)
# handles = [gzip.open(os.path.join(OUTPUT_Ecosystem_PATH, f"{cond}_vir.fasta.gz"), "wt") for cond in conds]
# for record in SeqIO.parse(gzip.open(IMGVR_FASTA_PATH, 'rt'), "fasta"):
#     UVIG = record.id.split("|")[0]
#     for i, cond in enumerate(conds_sets):
#         if UVIG in cond: SeqIO.write(record, handles[i], "fasta")

from Bio import SeqIO
import pandas as pd
import os


VHDB_FASTA_PATH = "/data/DeepMicroClass/data/virus-host-db/virushostdb.formatted.genomic.fna.gz"
DB_METADATA_PATH = "/data/DeepMicroClass/data/virus-host-db/virushostdb.tsv"

PROK_VIR_DIR = "/data/DeepMicroClass/data/virus-host-db/prok_vir"
os.makedirs(PROK_VIR_DIR, exist_ok=True)
EUK_VIR_DIR = "/data/DeepMicroClass/data/virus-host-db/euk_vir"
os.makedirs(EUK_VIR_DIR, exist_ok=True)

metadata = pd.read_table(DB_METADATA_PATH)
metadata.dropna(subset=["host tax id", "host lineage"], inplace=True)
metadata.reset_index(drop=True, inplace=True)

prok_vir = set()
euk_vir = set()
for i, row in metadata.iterrows():
    if row["host lineage"].split(";")[0] in {"Archaea", "Bacteria"}:
        for virus in row["refseq id"].split(", "): prok_vir.add(virus)
    elif row["host lineage"].split(";")[0] == "Eukaryota":
        for virus in row["refseq id"].split(", "): euk_vir.add(virus)

for record in SeqIO.parse(VHDB_FASTA_PATH, "fasta"):
    if record.id in prok_vir:
        SeqIO.write(record, os.path.join(PROK_VIR_DIR, f"{record.id}.fasta"), "fasta")
    elif record.id in euk_vir:
        SeqIO.write(record, os.path.join(EUK_VIR_DIR, f"{record.id}.fasta"), "fasta")


# Archaea_VIR_DIR = "/data/DeepMicroClass/data/virus-host-db/Archaea_vir"
# os.makedirs(Archaea_VIR_DIR, exist_ok=True)
# Bacteria_VIR_DIR = "/data/DeepMicroClass/data/virus-host-db/Bacteria_vir"
# os.makedirs(Bacteria_VIR_DIR, exist_ok=True)

# Archaea_vir = set()
# Bacteria_vir = set()
# for i, row in metadata.iterrows():
#     if row["host lineage"].split(";")[0] == "Archaea":
#         # Archaea_vir.update(row["refseq id"].split(", "))
#         for virus in row["refseq id"].split(", "): Archaea_vir.add(virus)
#     elif row["host lineage"].split(";")[0] == "Bacteria":
#         # Bacteria_vir.update(row["refseq id"].split(", "))
#         for virus in row["refseq id"].split(", "): Bacteria_vir.add(virus)

# for record in SeqIO.parse(VHDB_FASTA_PATH, "fasta"):
#     if record.id in Archaea_vir:
#         SeqIO.write(record, os.path.join(Archaea_VIR_DIR, f"{record.id}.fasta"), "fasta")
#     elif record.id in Bacteria_vir:
#         SeqIO.write(record, os.path.join(Bacteria_VIR_DIR, f"{record.id}.fasta"), "fasta")

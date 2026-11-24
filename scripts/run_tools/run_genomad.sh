#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate genomad
source ~/Downloads/mmseqs/util/bash-completion.sh

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    # python3 genomad/cli.py annotate --cleanup --sensitivity 4.2 --evalue 0.001 --threads 24 --splits 8 --lenient-taxonomy --full-ictv-lineage $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad /data/genomad_db
    # python3 genomad/cli.py marker-classification --threads 24 $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad /data/genomad_db
    # python3 genomad/cli.py nn-classification --cleanup --threads 24 $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad
    # python3 genomad/cli.py aggregated-classification $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad
    # python3 genomad/cli.py score-calibration $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad
    # python3 genomad/cli.py summary $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad
done

# micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
# python3 scripts/parse_result/parse_genomad.py


# --------------------- #
# CRC contig prediction #
# --------------------- #
for f in /storage1/liucongcong/crc/{assembly,fastq/OEP001340}/*.assembly.gz
do
    filename=$(basename -- $f)
    extension="${filename#*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python3 genomad/cli.py annotate --cleanup --sensitivity 4.2 --evalue 0.001 --threads 96 --splits 8 --lenient-taxonomy --full-ictv-lineage $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad /home/xuzhen/genomad_db
    python3 genomad/cli.py find-proviruses --cleanup --sensitivity 8.2 --crf-threshold 0.4 --marker-threshold 12.0 --evalue 0.001 --threads 96 --lenient-taxonomy --full-ictv-lineage $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad /home/xuzhen/genomad_db
    python3 genomad/cli.py marker-classification --threads 96 $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad /home/xuzhen/genomad_db
    python3 genomad/cli.py nn-classification --cleanup --threads 96 $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad
    python3 genomad/cli.py aggregated-classification $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad
    python3 genomad/cli.py score-calibration $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad
    python3 genomad/cli.py summary $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_geNomad
done

source ~/.bashrc && mamba activate MetaCAT
for f in /storage1/liucongcong/crc/{assembly,fastq/OEP001340}/*.assembly.gz
do
  filename=$(basename -- $f)
  extension="${filename##*.}"
  filename="${filename%%.*}"
  echo "Predicting $f"
  MetaCAT seed --threads 96 --fasta $f --output /storage1/xuzhen/DeepMetaClass_data/CRC_assembly_seed/${filename}.seed
  MetaCAT cluster --threads 96 --min-sequence-length 500 --fasta $f --coverage /storage1/xuzhen/DeepMetaClass_data/CRC_assembly_coverage/${filename}.coverage --output /storage1/xuzhen/DeepMetaClass_data/CRC_assembly_cluster_m500/${filename}.MetaCAT
done

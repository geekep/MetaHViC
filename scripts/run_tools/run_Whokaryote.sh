#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate Whokaryote
mkdir -p /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    whokaryote.py --contigs $f --outdir "/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote/${filename}" --model S --threads 24 --minsize 1
    # whokaryote.py --contigs $f --outdir "/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote/${filename}" --model T --threads 24 --minsize 1 --prodigal_file "/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote/${filename}/contigs_genes.gff"
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_Whokaryote.py

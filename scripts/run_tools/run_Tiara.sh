#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate Tiara

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    tiara -i $f -o "/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Tiara/${filename}.txt" --threads 24 --first_stage_kmer 6 --second_stage_kmer 7 --min_len 0 --probabilities --verbose # -p 0.2 0.3
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_Tiara.py

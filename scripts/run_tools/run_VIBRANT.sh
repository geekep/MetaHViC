#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate VIBRANT

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python3 /data/VIBRANT/VIBRANT_run.py -t 24 -f nucl -l 1000 -o 4 -i $f -folder /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_VIBRANT/
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_VIBRANT.py

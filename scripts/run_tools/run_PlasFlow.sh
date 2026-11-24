#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate PlasFlow
mkdir -p /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PlasFlow

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python3 /data/PlasFlow/PlasFlow.py --threshold 0.7 --input $f --output /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PlasFlow/$filename
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_PlasFlow.py

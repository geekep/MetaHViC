#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate DeepVirFinder

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python3 /data/DeepVirFinder/dvf.py -i $f -o /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_DeepVirFinder/ --core 12 --len 1
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_DeepVirFinder.py

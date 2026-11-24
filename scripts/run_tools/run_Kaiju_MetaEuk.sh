#!/bin/bash

eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate Kaiju && cd /media/xuzhen/Elements/DeepMetaClass

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    filename=$(basename -- "$f")
    extension="${filename#*.}"
    filename="${filename%%.*}"
    dir=$(dirname "$f")
    echo "Predicting $f"
    scripts/run_tools/Kaiju.sh --input_assembly $f --output_dir $dir --prefix $filename --threads 24
    scripts/run_tools/MetaEuk.sh --input_assembly $f --output_dir /home/xuzhen/Downloads/sequence_subsampled_1000_contig_6000_class_5 --prefix $filename --threads 24
done

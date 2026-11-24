#!/bin/bash

eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate Whokaryote && cd /media/xuzhen/Elements/DeepMetaClass

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    scripts/run_tools/parallel_prodigal.sh $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote/${filename} 24 /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote/${filename}/tmp
done

#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate VirSorter2

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    virsorter run -w /home/xuzhen/Downloads/sequence_subsampled_1000_contig_6000_class_5_Virsorter2/${filename} -i $f -j 24 --use-conda-off --provirus-off --max-orf-per-seq 10 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae all
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_VirSorter2.py

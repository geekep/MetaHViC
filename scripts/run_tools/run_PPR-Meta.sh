#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate PPR-Meta
export LD_LIBRARY_PATH=/usr/local/MATLAB/MATLAB_Runtime/v94/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v94/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v94/sys/os/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v94/extern/bin/glnxa64
mkdir -p /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PPR-Meta
cd /data/PPR-Meta

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    ./PPR_Meta $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PPR-Meta/${filename}.csv -t 0 -b 1000
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_PPR-Meta.py

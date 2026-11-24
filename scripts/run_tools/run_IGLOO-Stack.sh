#!/bin/bash

# Activate micromamba environment
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate genomad
source ~/Downloads/mmseqs/util/bash-completion.sh
mkdir -p /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_IGLOO-Stack_undersampling_nb_stacks_1_fungi

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_fungi/*.fa.gz
do
    dirpath=$(dirname -- "$f")
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    echo "Predicting $f"
    python3 genomad/cli.py host-classification --cleanup --threads 24 $f /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_IGLOO-Stack_undersampling_nb_stacks_1_fungi
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_IGLOO-Stack.py


# --------------------- #
# IMG contig prediction #
# --------------------- #
source ~/mmseqs/util/bash-completion.sh
mkdir -p /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack
mkdir -p /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack_prokaryotic
for f in /storage1/liucongcong/img/{assembly,img200/assembly-97}/img-*.assembly.gz
do
    dirpath=$(dirname -- "$f")
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%%.*}"
    if [ -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack/${filename}_host_classification/${filename}_host_classification.tsv" ]; then
        continue
    fi
    echo "Predicting $f"
    python3 genomad/cli.py host-classification --cleanup --threads 96 --batch-size 256 $f /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack
    python3 scripts/parse_result/extract_IGLOO-Stack.py extract --tsv /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack/${filename}_host_classification/${filename}_host_classification.tsv \
                                                                --fasta $f \
                                                                --class prokaryote \
                                                                --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_IGLOO_stack_prokaryotic/${filename}.prokaryote.assembly.fasta.gz
done


# --------------------- #
# CRC contig prediction #
# --------------------- #
source ~/mmseqs/util/bash-completion.sh
mkdir -p /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_IGLOO_stack
mkdir -p /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_IGLOO_stack_Prokaryote
for f in /storage1/liucongcong/crc/{assembly,fastq/OEP001340}/*.assembly.gz
do
    dirpath=$(dirname -- "$f")
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%%.*}"
    echo "Predicting $f"
    python3 genomad/cli.py host-classification --cleanup --threads 96 --batch-size 256 $f /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_IGLOO_stack
    python3 scripts/parse_result/extract_IGLOO-Stack.py extract --tsv /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_IGLOO_stack/${filename}_host_classification/${filename}_host_classification.tsv \
                                                                --fasta $f \
                                                                --class prokaryote \
                                                                --output /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_IGLOO_stack_Prokaryote/${filename}.prokaryote.assembly.fasta.gz
done

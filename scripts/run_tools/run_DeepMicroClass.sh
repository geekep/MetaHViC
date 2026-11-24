#!/bin/bash

eval "$(conda shell.bash hook)"
eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate DeepMicroClass

for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
do
    echo "Predicting $f"
    python3 DeepMicroClass/__main__.py predict --input_file $f --output_dir /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_DeepMicroClass --encoding onehot --mode hybrid --cpu_thread 12
    # python3 scripts/exp_script/DeepMicroClass_predict_pytorch.py --inputFile $f --outputDir /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_DeepMicroClass --modelDir DeepMicroClass/checkpoint --encoding onehot --mode hybrid
done

micromamba activate base && cd /media/xuzhen/Elements/DeepMetaClass
python3 scripts/parse_result/parse_DeepMicroClass.py


# ---------------------- #
# SPOT contig prediction #
# ---------------------- #
DeepMicroClass_RESULTS=/home/xuzhen/DeepMetaClass/data/SPOT_assembly_class_5 && mkdir -p ${DeepMicroClass_RESULTS}
model_dir="/home/xuzhen/DeepMetaClass/training_checkpoints"
for f in `ls /storage1/xuzhen/SPOT/assembly/*.fa`; do
    BASENAME=$(basename $f)
    FILESTEM="${BASENAME%.fa}"

    if [[ $FILESTEM =~ "SRZ1874" ]]; then
      if [[ ! -f ${DeepMicroClass_RESULTS}/${BASENAME}_pred_onehot_hybrid.tsv ]]; then
        python3 /home/xuzhen/DeepMicroClass/predict_pytorch.py --inputFile $f --outputDir ${DeepMicroClass_RESULTS} --modelDir ${model_dir} --encoding onehot --mode hybrid
        echo "${FILESTEM} is predicted successfully."
      else
        echo "${FILESTEM} has been predicted previously."
      fi
    fi

done


# --------------------- #
# IMG contig prediction #
# --------------------- #
for f in /storage1/liucongcong/img/{assembly,img200/assembly-97}/img-*.assembly.gz
do
    dirpath=$(dirname -- "$f")
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    if [ -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5/${filename}_pred_onehot_hybrid.tsv" ]; then
        continue
    fi
    echo "Predicting $f"
    python3 DeepMicroClass/__main__.py predict --input_file $f --output_dir /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5 --encoding onehot --mode hybrid --cpu_thread 32
    python3 DeepMicroClass/__main__.py extract --fasta $f --tsv /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5/${filename}_pred_onehot_hybrid.tsv --class Prokaryote,Plasmid,ProkaryoteVirus --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename%%.*}.prokaryotic.assembly.gz
    awk -F'\t' 'NR>1 && ($7!="Prokaryote" && $7!="Plasmid" && $7!="ProkaryoteVirus") {print $1}' /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5/${filename}_pred_onehot_hybrid.tsv > /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename%%.*}.non-prokaryotic.list
    awk -F'\t' 'NR>1 && ($7=="Eukaryote" || $7=="EukaryoteVirus") {print $1}' /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5/${filename}_pred_onehot_hybrid.tsv > /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename%%.*}.non-prokaryotic.list
done

parallel -j 64 '
filename=$(basename -- "{}")
[ -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename%%.*}.prokaryotic.assembly.gz" ] && exit 0
python3 DeepMicroClass/__main__.py extract \
  --fasta "{}" \
  --tsv /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5/${filename}_pred_onehot_hybrid.tsv \
  --class Prokaryote,Plasmid,ProkaryoteVirus \
  --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename%%.*}.prokaryotic.assembly.gz
' ::: /storage1/liucongcong/img/{assembly,img200/assembly-97}/img-*.assembly.gz

source /home/liucongcong/.bashrc
MetaCAT seed --threads 96 --fasta XXX.fa.gz --output XXX.seed
MetaCAT cluster --threads 96 --min-sequence-length 500 --fasta XXX.fa.gz --coverage XXX.coverage --seed XXX.seed --output XXX.MetaCAT
MetaCAT checkm2 --threads 96 --fasta-suffix fasta --fasta $(dirname -- "$XXX.001.fasta") --output XXX.checkm2
awk -i inplace -F'\t' -v OFS='\t' '{gsub(/\.MetaCAT\.cluster\./, ".MetaCAT-1.", $1)} 1' XXX.checkm2
tail -q -n +2 XXX1.checkm2 XXX2.checkm2 >> XXX3.checkm2
MetaCAT benchmarkRW --checkm2 XXX.checkm2 --output XXX --program MetaCAT MetaCAT-1

parallel -j 8 '
filename=$(basename -- "{}")
filename=${filename%%.*}
[ -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed" ] && exit 0
MetaCAT seed --threads 16 --fasta "{}" --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed
MetaCAT cluster --threads 16 --min-sequence-length 500 --fasta "{}" --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed_cluster_m500/${filename}.prokaryotic.MetaCAT
' ::: /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/img-*.prokaryotic.assembly.gz

ssh xuzhen@10.181.8.92 && mamba activate MetaCAT
for f in /storage1/liucongcong/img/{assembly,img200/assembly-97}/img-*.assembly.gz
do
  filename=$(basename -- $f)
  filename=${filename%%.*}
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed" ] && MetaCAT seed --threads 64 --fasta $f --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m500/${filename}.MetaCAT.mapping" ] && MetaCAT cluster --threads 64 --min-sequence-length 500 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m500/${filename}.MetaCAT
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m1000/${filename}.MetaCAT.mapping" ] && MetaCAT cluster --threads 64 --min-sequence-length 1000 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m1000/${filename}.MetaCAT
done
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 96 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m500 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 96 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_cluster_m1000 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2


for f in /storage1/liucongcong/img/{assembly,img200/assembly-97}/img-*.assembly.gz
do
  filename=$(basename -- $f)
  filename=${filename%%.*}
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500_NEW/${filename}.prokaryotic.MetaCAT.mapping" ] && MetaCAT cluster --threads 96 --min-sequence-length 500 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500_NEW/${filename}.prokaryotic.MetaCAT --exclude /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename}.non-prokaryotic.list
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000_NEW/${filename}.prokaryotic.MetaCAT.mapping" ] && MetaCAT cluster --threads 96 --min-sequence-length 1000 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_seed/${filename}.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000_NEW/${filename}.prokaryotic.MetaCAT --exclude /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/${filename}.non-prokaryotic.list
done
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 128 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500_NEW --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m500.checkm2
awk -i inplace -F'\t' -v OFS='\t' '{gsub(/\.prokaryotic\.MetaCAT\./, ".MetaCAT-2.", $1)} 1' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m500.checkm2
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 128 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000_NEW --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m1000.checkm2
awk -i inplace -F'\t' -v OFS='\t' '{gsub(/\.prokaryotic\.MetaCAT\./, ".MetaCAT-2.", $1)} 1' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m1000.checkm2

tail -q -n +2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m500.checkm2 >> /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2
MetaCAT benchmarkRW --checkm2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2_benchmarkRW/IMG.MetaCAT-m500 --program MetaCAT MetaCAT-1 MetaCAT-2
tail -q -n +2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-2-m1000.checkm2 >> /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2
MetaCAT benchmarkRW --checkm2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2_benchmarkRW/IMG.MetaCAT-m1000 --program MetaCAT MetaCAT-1 MetaCAT-2

sed -i '/MetaCAT-1/d' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2
sed -i '/MetaCAT-1/d' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2

for f in /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic/img-*.prokaryotic.assembly.gz
do
  filename=$(basename -- $f)
  filename=${filename%%.*}
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed" ] && MetaCAT seed --threads 96 --fasta $f --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500/${filename}.prokaryotic.MetaCAT.mapping" ] && MetaCAT cluster --threads 96 --min-sequence-length 500 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500/${filename}.prokaryotic.MetaCAT
  [ ! -f "/home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000/${filename}.prokaryotic.MetaCAT.mapping" ] && MetaCAT cluster --threads 96 --min-sequence-length 1000 --fasta $f --coverage /storage1/liucongcong/img/img200/img200.processed-files/img.metacat-1.0.3.processed-files/${filename}.metacat-1.0.3.coverage --seed /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_seed/${filename}.prokaryotic.seed --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000/${filename}.prokaryotic.MetaCAT
done
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 96 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m500 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m500.checkm2
awk -i inplace -F'\t' -v OFS='\t' '{gsub(/\.prokaryotic\.MetaCAT\./, ".MetaCAT-1.", $1)} 1' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m500.checkm2
MetaCAT checkm2 --checkm2 /home/liucongcong/program/checkm2-1.0.1/bin/checkm2 --threads 96 --fasta-suffix fasta --fasta /home/xuzhen/DeepMetaClass/data/IMG_assembly_class_5_prokaryotic_cluster_m1000 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m1000.checkm2
awk -i inplace -F'\t' -v OFS='\t' '{gsub(/\.prokaryotic\.MetaCAT\./, ".MetaCAT-1.", $1)} 1' /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m1000.checkm2

tail -q -n +2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m500.checkm2 >> /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2
MetaCAT benchmarkRW --checkm2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m500.checkm2 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2_benchmarkRW/IMG.MetaCAT-m500 --program MetaCAT MetaCAT-1
tail -q -n +2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-1-m1000.checkm2 >> /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2
MetaCAT benchmarkRW --checkm2 /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2/IMG.MetaCAT-m1000.checkm2 --output /home/xuzhen/DeepMetaClass/data/IMG_assembly_checkm2_benchmarkRW/IMG.MetaCAT-m1000 --program MetaCAT MetaCAT-1


# --------------------- #
# CRC contig prediction #
# --------------------- #
for f in /storage1/liucongcong/crc/{assembly,fastq/OEP001340}/*.assembly.gz
do
    dirpath=$(dirname -- "$f")
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    echo "Predicting $f"
    python3 DeepMicroClass/__main__.py predict --input_file $f --output_dir /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5 --encoding onehot --mode hybrid --cpu_thread 8
    python3 DeepMicroClass/__main__.py extract --fasta $f --tsv /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5/${filename}_pred_onehot_hybrid.tsv --class Prokaryote,Plasmid,ProkaryoteVirus --output /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5_prokaryotic/${filename%%.*}.prokaryotic.assembly.gz
done

parallel -j 64 '
filename=$(basename -- "{}")
python3 DeepMicroClass/__main__.py extract \
  --fasta "{}" \
  --tsv /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5/${filename}_pred_onehot_hybrid.tsv \
  --class Prokaryote,Plasmid,ProkaryoteVirus \
  --output /home/xuzhen/DeepMetaClass/data/CRC_assembly_class_5/${filename%%.*}.prokaryotic.assembly.gz
' ::: /storage1/liucongcong/crc/{assembly,fastq/OEP001340}/*.assembly.gz

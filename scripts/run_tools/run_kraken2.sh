#!/bin/bash

eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate kraken2 && cd /media/xuzhen/Elements/DeepMetaClass
export KRAKEN2_NUM_THREADS=24
export KRAKEN2_DB_PATH=/data/krakendb

# for f in /data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5/*.fa.gz
# do
#     filename=$(basename -- "$f")
#     extension="${filename#*.}"
#     filename="${filename%%.*}"
#     dir=$(dirname "$f")
#     echo "Predicting $f"
#     scripts/run_tools/kraken2.sh --input_assembly $f --output_dir $dir --prefix $filename --threads 24 --database k2_viral_20251015
# done

# These files can be found in $DBNAME/taxonomy/
# If not using custom sequences and not using one of the plasmid or non-redundant database libraries, skip downloading of the accession number to taxon maps
kraken2-build --download-taxonomy --skip-maps --db $DBNAME
kraken2-build --download-library bacteria --db $DBNAME
kraken2-build --download-library nr --protein --db $DBNAME
kraken2-build --add-to-library chr1.fa --db $DBNAME
# 35 and 31 for nucleotide database, 15 and 12 for protein database
# s = 7 for nucleotide database, s = 0 for protein database
# --no-masking option will skip masking of low-complexity sequences
kraken2-build --build --threads 24 --kmer-len 35 --minimizer-len 31 --minimizer-spaces 7 --no-masking --db $DBNAME
# remove intermediate files from the database directory
kraken2-build --clean
kraken2 --db $DBNAME --threads 56 --use-names --gzip-compressed --memory-mapping --report XXX.k2report --use-mpa-style --unclassified-out XXX.unclassified --classified-out XXX.classified --output XXX.kraken2 XXX.gz

# ------------------- #
# Microbiome Analysis #
# ------------------- #
working_dir=/data/m_samples && cd ${working_dir} && kraken2_RESULTS=kraken2_results && mkdir -p ${kraken2_RESULTS}
cat SraAccList_MicrobiomeAnalysis.txt | grep -v '^$' | while read id; do
    kraken2 --db k2_standard_eupath_20201202 --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/${id}.k2report --paired ${id}_1.fastq ${id}_2.fastq > ${kraken2_RESULTS}/${id}.kraken2
done

# cat SraAccList_MicrobiomeAnalysis.txt | grep -v '^$' | \
# parallel -j 1 "kraken2 --db k2_standard_eupath_20201202 --threads 24 \
#                         --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/{}.k2report \
#                         --paired {}_1.fastq {}_2.fastq > ${kraken2_RESULTS}/{}.kraken2"


# ----------------------- #
# Pathogen Identification #
# ----------------------- #
working_dir=/data/p_samples && cd ${working_dir} && kraken2_RESULTS=kraken2_results && mkdir -p ${kraken2_RESULTS}
export KRAKEN2_DEFAULT_DB=k2_standard_eupath_20201202
cat SraAccList_PathogenIdentification.txt | grep -v '^$' | while read id; do
    kraken2 --db k2_standard_eupath_20201202 --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/${id}.k2report --paired ${id}_1.fastq ${id}_2.fastq > ${kraken2_RESULTS}/${id}.kraken2
done

# cat SraAccList_PathogenIdentification.txt | grep -v '^$' | \
# parallel -j 1 "kraken2 --db k2_standard_eupath_20201202 --threads 24 \
#                        --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/{}.k2report \
#                        --paired {}_1.fastq {}_2.fastq > ${kraken2_RESULTS}/{}.kraken2"

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486971.kraken2 --include-children -s SRR12486971_1.fastq -s2 SRR12486971_2.fastq -t 723287 -r ${kraken2_RESULTS}/SRR12486971.k2report -o ${kraken2_RESULTS}/SRR12486971_A.algerae.tid723287.1.fa -o2 ${kraken2_RESULTS}/SRR12486971_A.algerae.tid723287.2.fa
bowtie2 -x /data/b_index/Aalgerae -f -p 8 -1 ${kraken2_RESULTS}/SRR12486971_A.algerae.tid723287.1.fa -2 ${kraken2_RESULTS}/SRR12486971_A.algerae.tid723287.2.fa -S ${kraken2_RESULTS}/SRR12486971_A.algerae_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486971_A.algerae_aligned.sam > ${kraken2_RESULTS}/SRR12486971_A.algerae_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486971_A.algerae_aligned.bam -o ${kraken2_RESULTS}/SRR12486971_A.algerae_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486971_A.algerae_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486972.kraken2 --include-children -s SRR12486972_1.fastq -s2 SRR12486972_2.fastq -t 5059 -r ${kraken2_RESULTS}/SRR12486972.k2report -o ${kraken2_RESULTS}/SRR12486972_A.flavus.tid5059.1.fa -o2 ${kraken2_RESULTS}/SRR12486972_A.flavus.tid5059.2.fa
bowtie2 -x /data/b_index/Aflavus -f -p 8 -1 ${kraken2_RESULTS}/SRR12486972_A.flavus.tid5059.1.fa -2 ${kraken2_RESULTS}/SRR12486972_A.flavus.tid5059.2.fa -S ${kraken2_RESULTS}/SRR12486972_A.flavus_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486972_A.flavus_aligned.sam > ${kraken2_RESULTS}/SRR12486972_A.flavus_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486972_A.flavus_aligned.bam -o ${kraken2_RESULTS}/SRR12486972_A.flavus_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486972_A.flavus_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486974.kraken2 --include-children -s SRR12486974_1.fastq -s2 SRR12486974_2.fastq -t 5476 -r ${kraken2_RESULTS}/SRR12486974.k2report -o ${kraken2_RESULTS}/SRR12486974_C.albicans.tid5476.1.fa -o2 ${kraken2_RESULTS}/SRR12486974_C.albicans.tid5476.2.fa
bowtie2 -x /data/b_index/Calbicans -f -p 8 -1 ${kraken2_RESULTS}/SRR12486974_C.albicans.tid5476.1.fa -2 ${kraken2_RESULTS}/SRR12486974_C.albicans.tid5476.2.fa -S ${kraken2_RESULTS}/SRR12486974_C.albicans_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486974_C.albicans_aligned.sam > ${kraken2_RESULTS}/SRR12486974_C.albicans_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486974_C.albicans_aligned.bam -o ${kraken2_RESULTS}/SRR12486974_C.albicans_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486974_C.albicans_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486978.kraken2 --include-children -s SRR12486978_1.fastq -s2 SRR12486978_2.fastq -t 1774 -r ${kraken2_RESULTS}/SRR12486978.k2report -o ${kraken2_RESULTS}/SRR12486978_M.chelonae.tid1774.1.fa -o2 ${kraken2_RESULTS}/SRR12486978_M.chelonae.tid1774.2.fa
bowtie2 -x /data/b_index/Mchelonae -f -p 8 -1 ${kraken2_RESULTS}/SRR12486978_M.chelonae.tid1774.1.fa -2 ${kraken2_RESULTS}/SRR12486978_M.chelonae.tid1774.2.fa -S ${kraken2_RESULTS}/SRR12486978_M.chelonae_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486978_M.chelonae_aligned.sam > ${kraken2_RESULTS}/SRR12486978_M.chelonae_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486978_M.chelonae_aligned.bam -o ${kraken2_RESULTS}/SRR12486978_M.chelonae_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486978_M.chelonae_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486983.kraken2 --include-children -s SRR12486983_1.fastq -s2 SRR12486983_2.fastq -t 10298 -r ${kraken2_RESULTS}/SRR12486983.k2report -o ${kraken2_RESULTS}/SRR12486983_HSV1.tid10298.1.fa -o2 ${kraken2_RESULTS}/SRR12486983_HSV1.tid10298.2.fa
bowtie2 -x /data/b_index/HSV1 -f -p 8 -1 ${kraken2_RESULTS}/SRR12486983_HSV1.tid10298.1.fa -2 ${kraken2_RESULTS}/SRR12486983_HSV1.tid10298.2.fa -S ${kraken2_RESULTS}/SRR12486983_HSV1_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486983_HSV1_aligned.sam > ${kraken2_RESULTS}/SRR12486983_HSV1_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486983_HSV1_aligned.bam -o ${kraken2_RESULTS}/SRR12486983_HSV1_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486983_HSV1_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486988.kraken2 --include-children -s SRR12486988_1.fastq -s2 SRR12486988_2.fastq -t 61605 -r ${kraken2_RESULTS}/SRR12486988.k2report -o ${kraken2_RESULTS}/SRR12486988_A.lugunensis.tid61605.1.fa -o2 ${kraken2_RESULTS}/SRR12486988_A.lugunensis.tid61605.2.fa
bowtie2 -x /data/b_index/Alug -f -p 8 -1 ${kraken2_RESULTS}/SRR12486988_A.lugunensis.tid61605.1.fa -2 ${kraken2_RESULTS}/SRR12486988_A.lugunensis.tid61605.2.fa -S ${kraken2_RESULTS}/SRR12486988_A.lugunensis_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486988_A.lugunensis_aligned.sam > ${kraken2_RESULTS}/SRR12486988_A.lugunensis_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486988_A.lugunensis_aligned.bam -o ${kraken2_RESULTS}/SRR12486988_A.lugunensis_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486988_A.lugunensis_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486989.kraken2 --include-children -s SRR12486989_1.fastq -s2 SRR12486989_2.fastq -t 1311 -r ${kraken2_RESULTS}/SRR12486989.k2report -o ${kraken2_RESULTS}/SRR12486989_S.agalactiae.tid1311.1.fa -o2 ${kraken2_RESULTS}/SRR12486989_S.agalactiae.tid1311.2.fa
bowtie2 -x /data/b_index/Sagalactiae -f -p 8 -1 ${kraken2_RESULTS}/SRR12486989_S.agalactiae.tid1311.1.fa -2 ${kraken2_RESULTS}/SRR12486989_S.agalactiae.tid1311.2.fa -S ${kraken2_RESULTS}/SRR12486989_S.agalactiae_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486989_S.agalactiae_aligned.sam > ${kraken2_RESULTS}/SRR12486989_S.agalactiae_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486989_S.agalactiae_aligned.bam -o ${kraken2_RESULTS}/SRR12486989_S.agalactiae_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486989_S.agalactiae_sorted.bam

python /data/KrakenTools/extract_kraken_reads.py -k ${kraken2_RESULTS}/SRR12486990.kraken2 --include-children -s SRR12486990_1.fastq -s2 SRR12486990_2.fastq -t 1280 -r ${kraken2_RESULTS}/SRR12486990.k2report -o ${kraken2_RESULTS}/SRR12486990_S.aureus.tid1280.1.fa -o2 ${kraken2_RESULTS}/SRR12486990_S.aureus.tid1280.2.fa
bowtie2 -x /data/b_index/Saureus -f -p 8 -1 ${kraken2_RESULTS}/SRR12486990_S.aureus.tid1280.1.fa -2 ${kraken2_RESULTS}/SRR12486990_S.aureus.tid1280.2.fa -S ${kraken2_RESULTS}/SRR12486990_S.aureus_aligned.sam
samtools view -bS -F4 ${kraken2_RESULTS}/SRR12486990_S.aureus_aligned.sam > ${kraken2_RESULTS}/SRR12486990_S.aureus_aligned.bam
samtools sort ${kraken2_RESULTS}/SRR12486990_S.aureus_aligned.bam -o ${kraken2_RESULTS}/SRR12486990_S.aureus_sorted.bam
samtools index ${kraken2_RESULTS}/SRR12486990_S.aureus_sorted.bam


# ----------------------------------- #
# CAMI II Clinical pathogen detection #
# ----------------------------------- #
working_dir=/data/patmg_CAMI2 && cd ${working_dir} && kraken2_RESULTS=kraken2_results && mkdir -p ${kraken2_RESULTS}
kraken2 --db k2_eupathdb48_20201113 --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/patmg_CAMI2_short_read.k2report --paired --gzip-compressed patmg_CAMI2_short_read_R1.fastq.gz patmg_CAMI2_short_read_R2.fastq.gz > ${kraken2_RESULTS}/patmg_CAMI2_short_read.kraken2
kraken2 --db k2_eupathdb48_20230407 --minimum-hit-groups 3 --report-minimizer-data --report ${kraken2_RESULTS}/patmg_CAMI2_short_read.k2report --paired --gzip-compressed patmg_CAMI2_short_read_R1.fastq.gz patmg_CAMI2_short_read_R2.fastq.gz > ${kraken2_RESULTS}/patmg_CAMI2_short_read.kraken2

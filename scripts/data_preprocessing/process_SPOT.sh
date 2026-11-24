#!/bin/bash

eval "$(~/.local/bin/micromamba shell hook --shell bash)"
micromamba activate base

# ------------------- #
# mapping to assemlby #
# ------------------- #
read_dir="/data/ESP_AE_MG_clean_reads"
assembly=/data/assembly/ESP_AE_MG_AM_min_2000_newbler_toAmos_minimus2_id0.98_renamed.fa
MAPPING_RESULTS=/data/SPOT_results/02_mapping_results && mkdir -p ${MAPPING_RESULTS}

# align reads to assemblies with bwa default
threads=64
for fwd_read in `ls ${read_dir}/*_fwd.fq.gz`; do
    base_name=`basename ${fwd_read}`
    rev_read=${fwd_read/_fwd/_rev}
    output_prefix=${base_name%_fwd.fq.gz}
    echo ${output_prefix}
    # index
    if [[ ! -f ${assembly}.bwt ]]; then
        bwa index $assembly
    fi
    # align, view to bam and sort
    bwa mem -t $threads $assembly $fwd_read $rev_read \
            | samtools view -@ $threads -bS - \
            | samtools sort -T ./tmp-samtools -@ $threads -O bam -o ${MAPPING_RESULTS}/${output_prefix}.bam -
done

# filter alignment at 70% identity
threads=24
identity=0.7
for ALIGNED_BAM in `ls ${MAPPING_RESULTS}/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    bam_file=${ALIGNED_BAM}
    bam_file2=${MAPPING_RESULTS}/${BAM_PREFIX}_id_${identity}.bam
    RUST_BACKTRACE=1 /home/xuzhen/.cargo/bin/coverm filter --bam-files ${bam_file} --threads ${threads} --output-bam-files ${bam_file2} --min-read-percent-identity ${identity} 
done

# align reads to assemblies with segment at 70% identity
threads=64
identity=70
for fwd_read in `ls ${read_dir}/*_fwd.fq.gz`; do
    base_name=`basename ${fwd_read}`
    rev_read=${fwd_read/_fwd/_rev}
    output_prefix=${base_name%_fwd.fq.gz}
    echo ${output_prefix}
    # index
    if [[ ! -f ${assembly}.idx ]]; then
        segemehl.x -x ${assembly}.idx -d $assembly -t ${threads}
    fi
    # align, view to bam and sort
    segemehl.x --silent --database $assembly -i ${assembly}.idx --query ${fwd_read} --mate ${rev_read} --threads ${threads} --accuracy ${identity} --hitstrategy 1 \
	| samtools view -@ ${threads} -bS - \
	| samtools sort -T ./tmp-samtools -@ $threads -O bam -o ${MAPPING_RESULTS}/${output_prefix}.bam -
done


# --------- #
# MG bamcov #
# --------- #
READ_COUNT_RESULTS="/data/SPOT_results/03_bamcov_readcnts_ID70" && mkdir -p $READ_COUNT_RESULTS
for ALIGNED_BAM in `ls $MAPPING_RESULTS/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    # bamcov tsv format: rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    bamcov --output ${READ_COUNT_RESULTS}/${BAM_PREFIX}_bamcov.tsv --min-read-len 30 --min-MQ 0 --min-BQ 0 $ALIGNED_BAM 
    awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $3, $4}' ${READ_COUNT_RESULTS}/${BAM_PREFIX}_bamcov.tsv > ${READ_COUNT_RESULTS}/${BAM_PREFIX}_name2len2reads.tsv
done

# --------- #
# MG CoverM #
# --------- #
MG_coverM_RESULTS=/data/VirEncAct/ESPAE/03_coverM_readcnts_ID98 && mkdir -p ${MG_coverM_RESULTS}
threads=24
for ALIGNED_BAM in `ls $MAPPING_RESULTS/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    bam_file=${ALIGNED_BAM}
    count_file=${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM.tsv
    # convert bam to coverM count table
    # --methods variance cause some problem
    RUST_BACKTRACE=1 /home/xuzhen/.cargo/bin/coverm contig \
        --bam-files ${bam_file} \
        --threads ${threads} \
        --min-read-aligned-length 30 \
        --min-read-percent-identity 0.98 \
        --min-read-aligned-percent 0.90 \
        --min-read-aligned-length-pair 30 \
        --min-read-percent-identity-pair 0.98 \
        --min-read-aligned-percent-pair 0.90 \
        --proper-pairs-only \
        --output-format dense \
        --methods count length covered_bases covered_fraction reads_per_base mean trimmed_mean rpkm \
        --min-covered-fraction 0 \
        --contig-end-exclusion 75 \
        --no-zeros \
        --trim-min 0.01 \
        --trim-max 0.99 1> ${count_file} 2> ${count_file}.log
    # get name2len2reads
    awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $3, $2}' ${count_file} > ${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM_name2len2reads.tsv 
done

MG_coverM_RESULTS=/data/SPOT_results/03_coverM_readcnts_ID70_pair && mkdir -p ${MG_coverM_RESULTS}
MAPPING_RESULTS=/data/VirEncAct/ESPAE/02_mapping_results/bwa
for ALIGNED_BAM in `ls $MAPPING_RESULTS/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    bam_file=${ALIGNED_BAM}
    count_file=${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM.tsv
    # convert bam to coverM count table
    # --methods variance cause some problem
    echo "counting ${bam_file} ..."
    RUST_BACKTRACE=1 /home/shengwei/.cargo/bin/coverm contig \
        --bam-files ${bam_file} \
        --threads ${threads} \
        --min-read-aligned-length 30 \
        --min-read-percent-identity 0.70 \
        --min-read-aligned-percent 0 \
        --min-read-percent-identity-pair 0.70 \
        --min-read-aligned-percent-pair 0.90 \
        --proper-pairs-only \
        --output-format dense \
        --methods count length covered_bases covered_fraction reads_per_base mean trimmed_mean rpkm \
        --min-covered-fraction 0 \
        --contig-end-exclusion 75 \
        --trim-min 0.05 \
        --trim-max 0.95 1> ${count_file} 2> ${count_file}.log
    # get name2len2reads
    awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $3, $2}' ${count_file} > ${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM_name2len2reads.tsv 
done

MG_coverM_RESULTS=/data/SPOT_results/03_coverM_readcnts_ID70 && mkdir -p ${MG_coverM_RESULTS}
MAPPING_RESULTS=/data/SPOT_results/02_mapping_results/bwa
for ALIGNED_BAM in `ls $MAPPING_RESULTS/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    bam_file=${ALIGNED_BAM}
    count_file=${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM.tsv
    # convert bam to coverM count table
    # --methods variance cause some problem
    echo "counting ${bam_file} ..."
    RUST_BACKTRACE=1 /home/xuzhen/.cargo/bin/coverm contig \
        --bam-files ${bam_file} \
        --threads ${threads} \
        --min-read-aligned-length 30 \
        --min-read-percent-identity 0.70 \
        --min-read-aligned-percent 0 \
        --output-format dense \
        --methods count length covered_bases covered_fraction reads_per_base mean trimmed_mean rpkm \
        --min-covered-fraction 0 \
        --contig-end-exclusion 75 \
        --trim-min 0.05 \
        --trim-max 0.95 1> ${count_file} 2> ${count_file}.log
    # get name2len2reads
    awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $3, $2}' ${count_file} > ${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM_name2len2reads.tsv 
done


# ----------- #
# count reads #
# ----------- #
DeepMicroClass_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID98 && mkdir -p ${DeepMicroClass_COUNTS}
DeepMicroClass_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID70 && mkdir -p ${DeepMicroClass_COUNTS}

DeepMicroClass_RESULTS=/data/SPOT_results/01_DeepMicroClass
DeepMicroClass_Prediction_file=${DeepMicroClass_RESULTS}/"ESP_AE_MG_AM_min_2000_newbler_toAmos_minimus2_id0.98_renamed.fa_pred_one-hot_hybrid.txt"
DeepMicroClass_Prediction_file=${DeepMicroClass_RESULTS}/"ESP_AEAM_DuraAMPM_MG_min_2000_newbler_toAmos_minimus2_id0.98_renamed.fa_pred_one-hot_hybrid.txt"
DeepMicroClass_Prediction_file=${DeepMicroClass_RESULTS}/"ESP_AEAM_DuraAMPM_MG_min_2000_newbler_toAmos_minimus2_id0.98_renamed.fa_pred_one-hot_hybrid_contig2type.tsv"

for ALIGNED_BAM in `ls $MAPPING_RESULTS/*.bam`; do
    BAM_BASENAME=$(basename "$ALIGNED_BAM")
    BAM_PREFIX="${BAM_BASENAME%.bam}"
    bam_file=${ALIGNED_BAM}
    name_len_cnt_file=${MG_coverM_RESULTS}/${BAM_PREFIX}_coverM_name2len2reads.tsv
    echo -e "current sample is ${BAM_PREFIX}"
    for l in 500 1000 2000 3000 4000 5000; do
        echo -e "\t => current length is $l"
        /data/VirEncAct/ESPAE/scripts/calculate_SeqType_counts.py ${DeepMicroClass_Prediction_file} ${name_len_cnt_file} -m $l --output_dir ${DeepMicroClass_COUNTS} --prefix ${BAM_PREFIX} --force
    done
done

# ----------- #
# count reads #
# ----------- #
READ_COUNT_RESULTS=/data/VirEncAct/ESPAE/03_coverM_readcnts_ID98
READ_COUNT_RESULTS=/data/VirEncAct/ESPAE/03_coverM_readcnts_ID70

SeqType_COUNTS=/data/SPOT_results/04_SeqType_Counts_ID70_coverM && mkdir -p ${SeqType_COUNTS}
SeqType_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID98 && mkdir -p ${SeqType_COUNTS}
SeqType_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID70 && mkdir -p ${SeqType_COUNTS}

kaiju_contig2type="DeepMicroClass_kaiju_metaeuk/ESPDTS_kaiju.names_contig2type.tsv"
metaeuk_contig2type="DeepMicroClass_kaiju_metaeuk/metaeuk_tax_results_tax_per_contig_contig2type.tsv"
DeepMicroClass_contig2type=${DeepMicroClass_RESULTS}/"ESP_AEAM_DuraAMPM_MG_min_2000_newbler_toAmos_minimus2_id0.98_renamed.fa_pred_one-hot_hybrid_contig2type.tsv"
for name_len_cnt_file in `ls ${READ_COUNT_RESULTS}/*_name2len2reads.tsv`; do
    BASENAME=$(basename "$name_len_cnt_file")
    PREFIX="${BASENAME%_name2len2reads.tsv}"
    echo -e "current sample is ${PREFIX}"

    # for kaiju
    prefix=kaiju_${PREFIX}
    for l in 500 1000 2000 3000 4000 5000; do
        echo -e "\t => current length is $l"
        python3 ./scripts/calculate_SeqType_counts.py ${kaiju_contig2type} ${name_len_cnt_file} -m $l --output_dir ${SeqType_COUNTS} --prefix ${prefix} --force
    done
    
    # for metaeuk
    prefix=metaeuk_${PREFIX}
    for l in 500 1000 2000 3000 4000 5000; do
        echo -e "\t => current length is $l"
        python3 ./scripts/calculate_SeqType_counts.py ${metaeuk_contig2type} ${name_len_cnt_file} -m $l --output_dir ${SeqType_COUNTS} --prefix ${prefix} --force
    done

    # for DeepMicroClass
    prefix=DeepMicroClass_${PREFIX}
    for l in 500 1000 2000 3000 4000 5000; do
        echo -e "\t => current length is $l"
        python3 /data/VirEncAct/ESPAE/scripts/calculate_SeqType_counts.py ${DeepMicroClass_contig2type} ${name_len_cnt_file} -m $l --output_dir ${SeqType_COUNTS} --prefix ${prefix} --force
    done
done


# ------------ #
# merge counts #
# ------------ #
SeqType_COUNTS=/data/SPOT_results/04_SeqType_Counts_ID70_coverM
SeqType_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID98
SeqType_COUNTS=/data/SPOT_results/04_DeepMicroClass_Counts_ID70

MERGED_SeqType_COUNTS=/data/SPOT_results/05_MERGED_SeqType_Counts_ID98 && mkdir -p ${MERGED_SeqType_COUNTS}
MERGED_SeqType_COUNTS=/data/SPOT_results/05_MERGED_SeqType_Counts_ID70 && mkdir -p ${MERGED_SeqType_COUNTS}

metadata="plot_figures/plot_SPOT/ESP_AE_metadata_CTD.tsv"
for l in 500 1000 2000 3000 4000 5000; do
    for count_prefix in "DeepMicroClass_"; do
        python3 /home/xuzhen/VirEncAct/ESPAE/scripts/merge_SeqType_counts.py ${SeqType_COUNTS} ${count_prefix} \
		        --suffix _min_${l}_ReadCnt_per_SeqType.tsv \
		        --prefix ${count_prefix}_MERGED_min_${l} \
		        -o ${MERGED_SeqType_COUNTS}
    done
    #python3 ~/scripts/add_metadata_to_DeepMicroClass_pct_table.py ${MERGED_SeqType_COUNTS}/${count_prefix}_MERGED_min_${l}_
    #python3 ~/scripts/merge_DeepMicroClass_class_counts.py ${DeepMicroClass_COUNTS} --suffix _min_${l}_ContigCnt.tsv --prefix DeepMicroClass_MERGED_min_${l} -o ${MERGED_DeepMicroClass_COUNTS}
    #python3 ~/scripts/add_metadata_to_DeepMicroClass_pct_table.py ${MERGED_DeepMicroClass_COUNTS}/DeepMicroClass_MERGED_min_${l}_pct.tsv ${metadata} -o ${MERGED_DeepMicroClass_COUNTS} 
    #python3 ~/scripts/add_metadata_to_DeepMicroClass_pct_table.py ${MERGED_DeepMicroClass_COUNTS}/DeepMicroClass_MERGED_min_${l}_cnt.tsv ${metadata} -o ${MERGED_DeepMicroClass_COUNTS} 
    #python3 ~/scripts/add_metadata_to_DeepMicroClass_pct_table.py ${MERGED_DeepMicroClass_COUNTS}/DeepMicroClass_MERGED_min_${l}_pct.tsv.ContigCnt ${metadata} -o ${MERGED_DeepMicroClass_COUNTS} -p DeepMicroClass_MERGED_min_${l}_ContigCntPct -f
done

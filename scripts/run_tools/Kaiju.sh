#!/bin/bash

# --------- #
# FUNCTIONS #
# --------- #

# help
usage()
{
    echo -e "\nusage: $0 -i <input_assembly> -o <output_dir> -p <prefix> -t <threads>\n"
}

# Kaiju
kaiju_classify(){
    # arguments
    input_assembly=$1
    KAIJU_RESULTS=$2
    prefix=$3
    threads=$4    

    # database
    nr=/data/kaiju_db/nr/kaiju_db_nr.fmi
    nr_euk=/data/kaiju_db/nr_euk/kaiju_db_nr_euk.fmi
    nodes_dmp=/data/kaiju_db/nr_euk/nodes.dmp
    names_dmp=/data/kaiju_db/nr_euk/names.dmp

    # search
    kaiju -z ${threads} \
        -t ${nodes_dmp} \
        -f ${nr_euk} \
        -e 5 \
        -E 0.05 \
        -x \
        -i ${input_assembly} \
        -o ${KAIJU_RESULTS}/${prefix}_kaiju.out

    # to krona
    kaiju2krona -t ${nodes_dmp} -n ${names_dmp} -i ${KAIJU_RESULTS}/${prefix}_kaiju.out -o ${KAIJU_RESULTS}/${prefix}_kaiju.out.krona
    ktImportText -o ${KAIJU_RESULTS}/${prefix}_kaiju.out.html ${KAIJU_RESULTS}/${prefix}_kaiju.out.krona

    # summarize classification
    for rank in phylum class order family genus species; do
        kaiju2table -p -t ${nodes_dmp} -n ${names_dmp} -r ${rank} -o ${KAIJU_RESULTS}/${prefix}_kaiju_summary_${rank}.tsv ${KAIJU_RESULTS}/${prefix}_kaiju.out
    done

    # addTaxonNames
    kaiju-addTaxonNames -p -t ${nodes_dmp} -n ${names_dmp} -i ${KAIJU_RESULTS}/${prefix}_kaiju.out -o ${KAIJU_RESULTS}/${prefix}_kaiju.names.out
}

# ---- #
# MAIN #
# ---- #

if [ "$#" != "8" ]; then
    echo -e "\nERROR: Not enough arguments provided!"
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    opt="$1"
    case "$opt" in
        -i|--input_assembly )
            shift
            input_assembly=$1
            ;;
        -o | --output_dir )
            shift
            output_dir=$1
            ;;
        -p | --prefix )
            shift
            prefix=$1
            ;;
        -t | --threads )
            shift
            threads=$1
            ;;
        -h | --help )
            usage
            exit
            ;;
        * )
            usage
            exit 1
            ;;
    esac
    shift
done

# run Kaiju
echo "The input assembly is: $input_assembly"
KAIJU_RESULTS=${output_dir}_Kaiju && mkdir -p ${KAIJU_RESULTS} && echo "The output directory is: $KAIJU_RESULTS"
echo "The prefix of output file is: $prefix"
echo "The number of threads is: $threads"

kaiju_classify $input_assembly $KAIJU_RESULTS $prefix $threads

#!/bin/bash

# --------- #
# FUNCTIONS #
# --------- #

# help
usage()
{
    echo -e "\nusage: $0 -i <input_assembly> -o <output_dir> -p <prefix> -t <threads>\n"
}

# MetaEuk
metaeuk_classify(){
    # arguments
    input_assembly=$1
    MetaEuk_RESULTS=$2
    prefix=$3
    threads=$4

    # database
    REF=/data/metaeuk_db/MMETSP_zenodo_3247846_uniclust90_2018_08_seed_valid_taxids
    metaeuk=/data/metaeuk/bin/metaeuk
    mmseqs=/home/xuzhen/Downloads/mmseqs/bin/mmseqs

    # create query db
    $mmseqs createdb ${input_assembly} ${MetaEuk_RESULTS}/queryContigDB

    # predict cds 
    $metaeuk easy-predict ${MetaEuk_RESULTS}/queryContigDB \
        ${REF} \
        ${MetaEuk_RESULTS}/${prefix}.fa \
        ${MetaEuk_RESULTS}/tempFolder \
        -s 4 --threads $threads

    # predict taxa
    $metaeuk taxtocontig ${MetaEuk_RESULTS}/queryContigDB \
        ${MetaEuk_RESULTS}/${prefix}.fa.fas \
        ${MetaEuk_RESULTS}/${prefix}.fa.headersMap.tsv \
        ${REF} \
        ${MetaEuk_RESULTS}/${prefix}_metaeuk_tax_results \
        ${MetaEuk_RESULTS}/tempFolder \
        --majority 0.5 --tax-lineage 1 --lca-mode 2 --threads $threads

    # remove query db
    mv ${MetaEuk_RESULTS}/queryContigDB.lookup ${MetaEuk_RESULTS}/${prefix}_backup_queryContigDB.lookup
    rm ${MetaEuk_RESULTS}/queryContigDB*
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

# run MetaEuk
echo "The input assembly is: $input_assembly"
MetaEuk_RESULTS=${output_dir}_MetaEuk && mkdir -p ${MetaEuk_RESULTS} && echo "The output directory is: $MetaEuk_RESULTS"
echo "The prefix of output file is: $prefix"
echo "The number of threads is: $threads"

metaeuk_classify $input_assembly $MetaEuk_RESULTS $prefix $threads

#!/bin/bash

# --------- #
# FUNCTIONS #
# --------- #

# help
usage()
{
    echo -e "\nusage: $0 -i <input_assembly> -o <output_dir> -p <prefix> -t <threads> -d <database>\n"
}

# kraken2
kraken2_classify() {
    # arguments
    input_assembly=$1
    kraken2_RESULTS=$2
    prefix=$3
    threads=$4
    # database
    krakendb=$5

    # predict taxa
    kraken2 --db $krakendb 
            --threads $threads \
            --use-names \
            --gzip-compressed \
            --memory-mapping \
            --report ${kraken2_RESULTS}/${prefix}.k2report \
            --use-mpa-style \
            --unclassified-out ${kraken2_RESULTS}/${prefix}.unclassified \
            --classified-out ${kraken2_RESULTS}/${prefix}.classified \
            --output ${kraken2_RESULTS}/${prefix}.kraken2 \
            ${input_assembly}
}

# ---- #
# MAIN #
# ---- #

if [ "$#" != "10" ]; then
    echo -e "\nERROR: Not enough arguments provided!"
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    opt="$1"
    case "$opt" in
        -i | --input_assembly )
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
        -d | --database )
            shift
            database=$1
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

# run kraken2
echo "The input assembly is: $input_assembly"
KRAKEN2_RESULTS=${output_dir}_kraken2 && mkdir -p ${KRAKEN2_RESULTS} && echo "The output directory is: $KRAKEN2_RESULTS"
echo "The prefix of output file is: $prefix"
echo "The number of threads is: $threads"
echo "The database is: $database"

kraken2_classify $input_assembly $KRAKEN2_RESULTS $prefix $threads $database

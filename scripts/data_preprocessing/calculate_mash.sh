#!/bin/bash

# Process input options
while getopts "d:nt:" opt; do
    case $opt in
        d) working_dir=$OPTARG;;
        t) threads=$OPTARG;;
        n) no_sketch=1;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac
done

working_dir=$working_dir && cd $working_dir
sketch_dir=/data/DeepMicroClass/data/mash/sketch
if [ ! -d "$sketch_dir" ]; then
    mkdir -p $sketch_dir
fi
if [ ! $no_sketch ]; then
    # mash sketch -p $threads -i /data/DeepMicroClass/data/clean_sequence/pre20_archaea.fa.gz -o $sketch_dir/pre20_archaea
    # mash sketch -p $threads -i /data/DeepMicroClass/data/clean_sequence/post20_archaea.fa.gz -o $sketch_dir/post20_archaea

    # mash sketch -p $threads -i /data/DeepMicroClass/data/clean_sequence/pre20_bacteria.fa.gz -o $sketch_dir/pre20_bacteria
    # mash sketch -p $threads -i /data/DeepMicroClass/data/clean_sequence/post20_bacteria.fa.gz -o $sketch_dir/post20_bacteria

    # mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/pre20_eukaryote.fa.gz -o $sketch_dir/pre20_eukaryote
    # mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/post20_eukaryote.fa.gz -o $sketch_dir/post20_eukaryote

    # mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/MMETSP_Train_Validate.fasta.gz -o $sketch_dir/MMETSP_Train_Validate
    # mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/MMETSP_Test.fasta.gz -o $sketch_dir/MMETSP_Test

    mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/pre20_fungi.fa.gz -o $sketch_dir/pre20_fungi
    mash sketch -p $threads -i /data/DeepMicroClass/data/raw_sequence/post20_fungi.fa.gz -o $sketch_dir/post20_fungi

    # mash sketch -p $threads -i /data/DeepMicroClass/data/plsdb/pre20_plasmid.fa.gz -o $sketch_dir/pre20_plasmid
    # mash sketch -p $threads -i /data/DeepMicroClass/data/plsdb/post20_plasmid.fa.gz -o $sketch_dir/post20_plasmid

    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/euk_vir_Train_Validate.fasta.gz -o $sketch_dir/euk_vir_Train_Validate
    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/euk_vir_Test.fasta.gz -o $sketch_dir/euk_vir_Test

    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/prok_vir_Train_Validate.fasta.gz -o $sketch_dir/prok_vir_Train_Validate
    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/prok_vir_Test.fasta.gz -o $sketch_dir/prok_vir_Test

    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/Archaea_vir_Train_Validate.fasta.gz -o $sketch_dir/Archaea_vir_Train_Validate
    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/Archaea_vir_Test.fasta.gz -o $sketch_dir/Archaea_vir_Test

    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/Bacteria_vir_Train_Validate.fasta.gz -o $sketch_dir/Bacteria_vir_Train_Validate
    # mash sketch -p $threads -i /data/DeepMicroClass/data/virus-host-db/Bacteria_vir_Test.fasta.gz -o $sketch_dir/Bacteria_vir_Test
fi

dist_dir=/data/DeepMicroClass/data/mash/dist
if [ ! -d "$dist_dir" ]; then
    mkdir -p $dist_dir
fi
# mash dist -p $threads $sketch_dir/pre20_archaea.msh $sketch_dir/post20_archaea.msh -t > $dist_dir/archaea.tsv
# mash dist -p $threads $sketch_dir/pre20_bacteria.msh $sketch_dir/post20_bacteria.msh -t > $dist_dir/bacteria.tsv
# mash dist -p $threads $sketch_dir/pre20_eukaryote.msh $sketch_dir/post20_eukaryote.msh -t > $dist_dir/eukaryote.tsv
# mash dist -p $threads $sketch_dir/MMETSP_Train_Validate.msh $sketch_dir/MMETSP_Test.msh -t > $dist_dir/MMETSP.tsv
mash dist -p $threads $sketch_dir/pre20_fungi.msh $sketch_dir/post20_fungi.msh -t > $dist_dir/fungi.tsv
# mash dist -p $threads $sketch_dir/pre20_plasmid.msh $sketch_dir/post20_plasmid.msh -t > $dist_dir/plasmid.tsv

#!/bin/bash

working_dir=data && mkdir -p $working_dir && cd $working_dir

if [ ! -f "genome_updater.sh" ]; then
    wget --quiet --show-progress https://raw.githubusercontent.com/pirovc/genome_updater/master/genome_updater.sh
    chmod +x genome_updater.sh
fi

##########################################################################################
###   Download bacteria, archaea, and eukaryota genomes from NCBI RefSeq and GenBank   ###
##########################################################################################
# A function named combine_fasta with parameters $input_dir and $output_file
combine_fasta () {
    input_dir=$1
    output_file=$2
    # macOS
    # gzcat $input_dir/*.fa.gz | gzip -c > $output_file
    # find $input_dir -type f -name "*.gz" -exec gzcat {} + | gzip -c > $output_file
    # Linux
    find $input_dir -type f -name "*.gz" -exec zcat {} + | gzip -c > $output_file
}

# Download archaea from RefSeq and GenBank, then combine all archaea genomes into a file
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/pre20_archaea" -T '2157' -E 20191230 -R 50 && combine_fasta "raw_sequence/pre20_archaea" "raw_sequence/pre20_archaea.fa.gz"
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/post20_archaea" -T '2157' -D 20200101 -R 50 && combine_fasta "raw_sequence/post20_archaea" "raw_sequence/post20_archaea.fa.gz"
# Download bacteria from RefSeq and GenBank, then combine all bacteria genomes into a file
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/pre20_bacteria" -T '2' -E 20191230 -R 50 && combine_fasta "raw_sequence/pre20_bacteria" "raw_sequence/pre20_bacteria.fa.gz"
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/post20_bacteria" -T '2' -D 20200101 -R 50 && combine_fasta "raw_sequence/post20_bacteria" "raw_sequence/post20_bacteria.fa.gz"
# Download eukaoryta genomes from Kaiju and PR2, then combine all eukaryota genomes into a file
# 554915: Eukaryota|Amoebozoa
# 554296: Eukaryota|Apusozoa
# 3027:   Eukaryota|Cryptophyceae
# 33682:  Eukaryota|Discoba|Euglenozoa
# 33634:  Eukaryota|Sar|Stramenopiles
# 33630:  Eukaryota|Sar|Alveolata
# 543769: Eukaryota|Sar|Rhizaria
# 5752:   Eukaryota|Discoba|Heterolobosea
# 2611341:Eukaryota|Metamonada
# 2763:   Eukaryota|Rhodophyta(red algae)
# 3041:   Eukaryota|Viridiplantae|Chlorophyta(green algae)
# 38254:  Eukaryota|Glaucocystophyceae
# 260810: Viruses|Riboviria|Orthornavirae|Negarnaviricota|Polyploviricotina|Insthoviricetes|Articulavirales|Orthomyxoviridae|Alphainfluenzavirus|Alphainfluenzavirus influenzae|Influenza A virus|H1N1 subtype
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/pre20_eukaryote" -T '554915,554296,3027,33682,33634,33630,543769,5752,2611341,2763,3041,38254,260810' -E 20191230 -R 50 && combine_fasta "raw_sequence/pre20_eukaryote" "raw_sequence/pre20_eukaryote.fa.gz"
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/post20_eukaryote" -T '554915,554296,3027,33682,33634,33630,543769,5752,2611341,2763,3041,38254,260810' -D 20200101 -R 50 && combine_fasta "raw_sequence/post20_eukaryote" "raw_sequence/post20_eukaryote.fa.gz"
# Download fungi genomes from RefSeq and GenBank, then combine all fungi genomes into a file
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/pre20_fungi" -g "fungi" -E 20191230 -R 50 && combine_fasta "raw_sequence/pre20_fungi" "raw_sequence/pre20_fungi.fa.gz"
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/post20_fungi" -g "fungi" -D 20200101 -R 50 && combine_fasta "raw_sequence/post20_fungi" "raw_sequence/post20_fungi.fa.gz"
# Download virus genomes from RefSeq and GenBank, then combine all virus genomes into a file
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/pre20_virus" -g "viral" -E 20191230 -R 50 && combine_fasta "raw_sequence/pre20_virus" "raw_sequence/pre20_virus.fa.gz"
./genome_updater.sh -M ncbi -d "genbank" -l "Complete Genome,Chromosome" -A 1 -f "genomic.fna.gz" -m -t 24 -o "raw_sequence/post20_virus" -g "viral" -D 20200101 -R 50 && combine_fasta "raw_sequence/post20_virus" "raw_sequence/post20_virus.fa.gz"

##################################
###   Download virus-host db   ###
##################################
vhdb_dir=$working_dir/virus-host-db && mkdir -p $vhdb_dir && cd $vhdb_dir
wget -q --show-progress https://www.genome.jp/ftp/db/virushostdb/virushostdb.formatted.genomic.fna.gz && gzip -d virushostdb.formatted.genomic.fna.gz
wget -q --show-progress https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
# Split vhdb according to prok and euk
python3 scripts/data_preprocessing/process_vhdb.py
# cd $working_dir && ln -sf $vhdb_dir vhdb

##########################
###   Download PLSDB   ###
##########################
plsdb_dir=$working_dir/plsdb && mkdir -p $plsdb_dir && cd $plsdb_dir
wget -q --show-progress https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plasmids_meta.tar.bz2 && tar -xjvf plasmids_meta.tar.bz2
wget -q --show-progress https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2 && bzip2 -d plsdb.fna.bz2
# Split plsdb according to date 2020-01-01
python3 scripts/data_preprocessing/process_plsdb.py

################################################################
###   Download MMETSP (https://zenodo.org/records/1212585)   ###
################################################################
mmetsp_dir=$working_dir/MMETSP && mkdir -p $mmetsp_dir && cd $mmetsp_dir
zenodo_get -R 32 1212585

##############################
###   Download IMG/VR v4   ###
##############################
IMGVR_dir=$working_dir/IMGVR && mkdir -p $IMGVR_dir && cd $IMGVR_dir
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=xuzhen3039@stu.xjtu.edu.cn' --data-urlencode 'password=bn6*@R9C-xEa7!8' -c cookies > /dev/null
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/IMG_VR/download/_JAMO/632a3d7d2de7c323533dac10/README-high_confidence.txt' -b cookies > README-high_confidence.txt
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/IMG_VR/download/_JAMO/632a3d7d2de7c323533dac0c/IMGVR_all_Sequence_information-high_confidence.tsv' -b cookies > IMGVR_all_Sequence_information.tsv
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/IMG_VR/download/_JAMO/632a3d7d2de7c323533dac0e/IMGVR_all_Host_information-high_confidence.tsv' -b cookies > IMGVR_all_Host_information.tsv
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/IMG_VR/download/_JAMO/632a3d7d2de7c323533dac0a/IMGVR_all_proteins-high_confidence.faa.gz' -b cookies > IMGVR_all_proteins.faa.gz
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/IMG_VR/download/_JAMO/632a3d7d2de7c323533dac08/IMGVR_all_nucleotides-high_confidence.fna.gz' -b cookies > IMGVR_all_nucleotides.fna.gz

#!/bin/bash

# python3 scripts/data_preprocessing/remove_plasmid.py --input_dir data/raw_sequence --output_dir data/clean_sequence

# dir=data/clean_sequence_filtered_fasta_0.1 && mkdir -p $dir && echo $dir
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/clean_sequence/pre20_archaea.fa.gz --test data/clean_sequence/post20_archaea.fa.gz --dist data/mash/dist/archaea.tsv --output $dir/archaea.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/clean_sequence/pre20_bacteria.fa.gz --test data/clean_sequence/post20_bacteria.fa.gz --dist data/mash/dist/bacteria.tsv --output $dir/bacteria.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/clean_sequence/pre20_eukaryote.fa.gz --test data/clean_sequence/post20_eukaryote.fa.gz --dist data/mash/dist/eukaryote.tsv --output $dir/eukaryote.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/plasmid/pre20_plasmid.fa.gz --test data/plasmid/post20_plasmid.fa.gz --dist data/mash/dist/plasmid.tsv --output $dir/plasmid.fa.gz
python3 scripts/data_preprocessing/remove_similar_seq.py --train data/virus-host-db/prok_vir_Train_Validate.fasta.gz --test data/virus-host-db/prok_vir_Test.fasta.gz --dist data/mash/dist/prokvir.tsv --output $dir/prokvir.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/virus-host-db/Archaea_vir_Train_Validate.fasta.gz --test data/virus-host-db/Archaea_vir_Test.fasta.gz --dist data/mash/dist/Archaea_vir.tsv --output $dir/Archaea_vir.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/virus-host-db/Bacteria_vir_Train_Validate.fasta.gz --test data/virus-host-db/Bacteria_vir_Test.fasta.gz --dist data/mash/dist/Bacteria_vir.tsv --output $dir/Bacteria_vir.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train data/virus-host-db/euk_vir_Train_Validate.fasta.gz --test data/virus-host-db/euk_vir_Test.fasta.gz --dist data/mash/dist/eukvir.tsv --output $dir/eukvir.fa.gz
# python3 scripts/data_preprocessing/remove_similar_seq.py --train /data/DeepMicroClass/data/raw_sequence/MMETSP_Train_Validate.fasta.gz --test /data/DeepMicroClass/data/raw_sequence/MMETSP_Test.fasta.gz --dist data/mash/dist/MMETSP.tsv --output $dir/MMETSP.fa.gz
# zcat $dir/eukaryote.fa.gz $dir/MMETSP.fa.gz > $dir/euk_MMETSP.fa.gz
# zcat $dir/archaea.fa.gz $dir/bacteria.fa.gz > $dir/prokaryote.fa.gz

# python3 scripts/data_preprocessing/split_fasta.py --input data/clean_sequence/pre20_archaea.fa.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/pre20_archaea.fa_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/clean_sequence/pre20_bacteria.fa.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/pre20_bacteria.fa_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/clean_sequence/pre20_prokaryote.fa.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/pre20_prokaryote.fa_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/clean_sequence/pre20_eukaryote.fa.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/pre20_eukaryote.fa_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/plasmid/pre20_plasmid.fa.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/pre20_plasmid.fa_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/virus-host-db/Archaea_vir_Train_Validate.fasta.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/Archaea_vir_Train_Validate.fasta_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/virus-host-db/Bacteria_vir_Train_Validate.fasta.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/Bacteria_vir_Train_Validate.fasta_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/virus-host-db/prok_vir_Train_Validate.fasta.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/prok_vir_Train_Validate.fasta_single_fasta
# python3 scripts/data_preprocessing/split_fasta.py --input data/virus-host-db/euk_vir_Train_Validate.fasta.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir data/clean_sequence/euk_vir_Train_Validate.fasta_single_fasta
python3 scripts/data_preprocessing/split_fasta.py --input /data/DeepMicroClass/data/raw_sequence/MMETSP_Train_Validate.fasta.gz && python3 scripts/data_preprocessing/fasta_to_onehot.py --input_dir /data/DeepMicroClass/data/raw_sequence/MMETSP_Train_Validate.fasta_single_fasta

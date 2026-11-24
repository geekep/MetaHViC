import argparse
import gzip
import pandas as pd
from Bio import SeqIO

def classify_sequences(tsv_file):
    # read the TSV file into a dataframe
    df = pd.read_csv(tsv_file, sep="\t")
    # function to find the class with the highest score
    def get_class(row):
        class_scores = row[['prokaryote', 'prokaryoteVirus', 'plasmid', 'eukaryote', 'eukaryoteVirus']]
        highest_class = class_scores.idxmax()
        return highest_class
    # apply the function to get the class for each sequence
    df["Predicted"] = df.apply(get_class, axis=1)
    return dict(zip(df["seq_name"], df["Predicted"]))

def extract_sequences_by_class(fasta_file, sequence_classes, desired_classes, output_fasta):
    desired_classes = [cls.strip() for cls in desired_classes.split(',') if cls.strip()]
    if not desired_classes: return
    if fasta_file.endswith((".fa", ".fasta", ".fna")) and not fasta_file.startswith("."):
        seq_records = SeqIO.parse(open(fasta_file, 'rt'), "fasta")
    elif fasta_file.endswith(".gz") and not fasta_file.startswith("."):
        seq_records = SeqIO.parse(gzip.open(fasta_file, 'rt'), "fasta")
    desired_records = (record for record in seq_records if record.id in sequence_classes and sequence_classes[record.id] in desired_classes)
    with gzip.open(output_fasta, "wt") as output_handle:
        SeqIO.write(desired_records, output_handle, "fasta")

def extract_wrapper(tsv_file, fasta_file, desired_classes, output_fasta):
    # classify sequences based on the TSV file
    sequence_classes = classify_sequences(tsv_file)
    # extract sequences of the desired class from the FASTA file
    extract_sequences_by_class(fasta_file, sequence_classes, desired_classes, output_fasta)
    print(f"Sequences of class {desired_classes} in {fasta_file} have been extracted to {output_fasta}.")


def main():
    parser = argparse.ArgumentParser(prog="IGLOO-Stack", description="IGLOO-Stack: A deep learning framework for classifying metagenomic sequences")
    subparsers = parser.add_subparsers()

    # Extract subcommand
    parser_extract = subparsers.add_parser("extract", help="Extract sequences of a particular class from a FASTA file")
    parser_extract.add_argument("--tsv", "-t", dest="tsv_file", help="Path to the TSV classification file", required=True)
    parser_extract.add_argument("--fasta", "-f", dest="fasta_file", help="Path to the input FASTA file", required=True)
    parser_extract.add_argument("--class", "-c", dest="desired_classes", help="Desired classes to extract sequences. It could be 'prokaryote', 'prokaryoteVirus', 'plasmid', 'eukaryote', 'eukaryoteVirus' or 'prokaryote,plasmid'", required=True)
    parser_extract.add_argument("--output", "-o", dest="output_fasta", help="Path to the output FASTA file", required=True)
    parser_extract.set_defaults(func=extract_wrapper)

    args = parser.parse_args()
    # Convert args to a dictionary
    args_dict = vars(args)
    # Remove 'func' from the dictionary
    func = args_dict.pop("func")
    # Call the function with the remaining arguments
    func(**args_dict)


if __name__ == "__main__":
    main()

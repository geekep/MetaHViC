# MetaHViC
MetaHViC: a Host-informed Virus Classifier for Metagenomic contigs

## Feature
MetaHViC's primary goal is to identify viruses and plasmids in sequencing data (isolates, metagenomes, and metatranscriptomes). Assigning contigs to biologically host-informed categories, including
- prokaryotic chromosomes
- prokaryotic viruses
- plasmids
- eukaryotic chromosomes
- eukaryotic viruses

### Installation
With Mamba, you will create an environment for MetaHViC and activate it before being able to use it.

```
# Create an environment for MetaHViC
mamba create -n MetaHViC -f requirements.txt -c conda-forge -c bioconda
# Activate the MetaHViC environment
mamba activate MetaHViC
```

### Executing MetaHViC
In this example, we will use an *Klebsiella pneumoniae* genome ([GCF_009025895.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009025895.1/)) as input. You can use any FASTA file containing nucleotide sequences as input.

The command to execute MetaHViC is structured like this:

```
MetaHViC host-classification [OPTIONS] INPUT OUTPUT DATABASE
```

The results will be written inside the `MetaHViC_output` directory.

In this example, the results of MetaHViC's analysis will be written to the `MetaHViC_output` directory, which will look like this:

```
MetaHViC_output
├── GCF_009025895.1_host_classification
╰── GCF_009025895.1_host_classification.log
```

## Citing MetaHViC
If you use MetaHViC in your work, please consider citing its manuscript:

> [**A Host-informed Virus Classifier for Metagenomic contigs**]()
>
>
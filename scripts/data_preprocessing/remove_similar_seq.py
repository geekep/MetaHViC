from Bio import SeqIO
import pandas as pd
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--train", "-i1", dest="train")
parser.add_argument("--test", "-i2", dest="test")
parser.add_argument("--output", "-o", dest="output")
parser.add_argument("--dist", "-d", dest="dist")
parser.add_argument("--threshold", "-t", dest="threshold", default=0.1, type=float)
args = parser.parse_args()

dist = pd.read_table(args.dist, index_col=0, delimiter='\t', engine='c', na_filter=False, low_memory=False)
with gzip.open(args.train, 'rt') as handle:
    train_id = [record.id for record in SeqIO.parse(handle, 'fasta') if record.id in dist.columns]
with gzip.open(args.test, 'rt') as handle:
    test_id = [record.id for record in SeqIO.parse(handle, 'fasta') if record.id in dist.index]
train_test_dist = dist.loc[test_id, train_id].T
train_test_overlap = set(train_test_dist.columns[(train_test_dist < args.threshold).sum(axis=0) != 0]) | set(train_test_dist.index[(train_test_dist < args.threshold).sum(axis=1) != 0])
test_id = set(id for id in test_id if id not in train_test_overlap)
with gzip.open(args.test, 'rt') as read_handle, gzip.open(args.output, 'wt') as write_handle:
    test_records = (record for record in SeqIO.parse(read_handle, 'fasta') if record.id in test_id)
    SeqIO.write(test_records, write_handle, 'fasta')

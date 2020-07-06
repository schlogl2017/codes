import os
import gzip
import glob
import csv
import pandas as pd
import itertools
from collections import defaultdict, Counter


print(os.getcwd())


def get_kmers(sequence, k):
    return (sequence[i:i+k] for i in range(len(sequence) - k + 1))


def is_header(line):
    return line[0] == '>'


def parse_multi_fasta_file_compressed_or_not(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        fasta_iter = (it[1] for it in itertools.groupby(f, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequences = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequences


def get_kmers_count_splitted_fasta_files(filename, k):
    count = Counter()
    for name, seq in parse_multi_fasta_file_compressed_or_not(filename):
        count.update(get_kmers(seq, k))
    return count


print(glob.glob('split_data/*.fasta'))
files = [filename for filename in glob.glob('split_data/*.fasta')]
print(files)


k = 7
counters_uniprot = [get_kmers_count_splitted_fasta_files(filename, k) for filename in files]
total_sprot = sum(counters_uniprot, Counter())

print('All done!')
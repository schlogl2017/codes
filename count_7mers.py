import os
import sys
import gzip
import glob
import itertools
from collections import Counter


if len(sys.argv) < 3:
    print('Usage: count_mers.py <filedir/files> <k (integer representing the mers lengths')
    sys.exit(1)

filename = sys.argv[1]
k = sys.argv[2]

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


files = [filename for filename in glob.glob(sys.argv[1])


for f in files:
    counters_uniprot = [get_kmers_count_splitted_fasta_files(filename, k) for filename in files]
    total_sprot = sum(counters_uniprot, Counter())
    
    with open('split_data/uniprot_sprot_7mers.csv', 'w') as fout:
    fout.write('Kmers' + ',' + 'Count' + '\n')
    for k, v in total_sprot.items():
        fout.write(k + ',' + str(v) + '\n')


print('All done!')

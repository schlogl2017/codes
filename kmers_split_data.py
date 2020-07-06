import os
import gzip
import glob
import itertools
from collections import Counter


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
        count.update(get_kmers(k, seq))
    return count


files = [filename for filename in glob.iglob('split_data/*.fasta')]


files_1 = files[:10]
print(len(files_1))
files_2 = files[10:20]
print(len(files_2))
files_3 = files[20:]
print(len(files_3))


k = 7
print('The length of the kmers is {}'.format(k))


counters_1 = [get_kmers_count_splitted_fasta_files(filename, k) for filename in files_1]
sprot1 = sum(counters_1, Counter())

counters_2 = [get_kmers_count_splitted_fasta_files(filename, k) for filename in files_2]
sprot2 = sum(counters_2, Counter())

counters_3 = [get_kmers_count_splitted_fasta_files(filename, k) for filename in files_3]
sprot3 = sum(counters_3, Counter())

uniprot_seven_mers = sum([sprot1, sprot2, sprot3], Counter())

with open('seven_mers.csv', 'w') as fout:
    for (k, v) in uniprot_seven_mers.items():
        fout.write(k + ',' + str(v) + '\n')


print('Script Done!')







































#!/usr/bin/env python
# coding: utf-8


import sys
import gzip
import time
import numpy as np
import itertools as it


start_time = time.process_time()


if len(sys.argv) != 3:
    print('Usage: permutation.py <input file> <output_file>')
    sys.exit(1)


filename = sys.argv[1]
output_file = sys.argv[2]


def is_header(line):
    """Check if the first line is a header line in a fasta file"""
    return line[0] == '>'


def parse_multi_fasta_file_compressed_or_not(filename):
    """Yields the name(ID) and a sequence as a tuple
    arg = a fasta file that can be compressed as a gzip file
    or a non compressed fasta file.    
    """
    print('Starting reading the fasta file' + '\n')
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        fasta_iter = (it[1] for it in it.groupby(f, is_header))
        for name in fasta_iter:
            name = next(name)
            sequences = ''.join(seq.strip() for seq in next(fasta_iter))
            yield name.strip()[1:], sequences
    print('Finishing reading the fasta file!')


with gzip.open(output_file +'.faa.gz', 'wt') as fout:
    for i, (name, seq) in enumerate(parse_multi_fasta_file_compressed_or_not(filename)):
        fout.write('>' + 'rand_seq_%d' % i + '\n')
        fout.write(''.join(np.random.permutation(list(seq))) +'\n')

        
print('Script finished!')
print('Time of execution: {} seconds'.format(time.process_time() - start_time))

sys.exit()


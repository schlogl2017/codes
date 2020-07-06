#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 16:39:20 2020

@author: paulo
"""
import os
import gzip
import glob
import itertools
from collections import Counter
from Bio import SeqIO



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


def count_fasta_files(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        return sum(g for g, _ in itertools.groupby(f, key=is_header))


def split_fasta_batch(iterator, file_size):
    loop = True
    while loop:
        size = []
        while len(size) < file_size:
            try:
                loop = iterator.__next__()
            except StopIteration:
                loop = None
            if loop is None:
                break
            size.append(loop)
        if size:
            yield size


filename = 'data/uniprot_sprot.fasta.gz'
with gzip.open(filename, 'rt') as fh:
    for i, batch in enumerate(split_fasta_batch(SeqIO.parse(fh, 'fasta'), 20000), start=1):
        filename = 'split_data/uniprot_sprot_{}.fasta'.format(i)
        count = SeqIO.write(batch, filename, 'fasta')
        print('Wrote {} records to {}'.format(count, filename))

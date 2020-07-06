#!usr/bin/env python

import gzip
import itertools


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


def fasta_item_counter(filename):
    """It opens and check the number of items in the fasta file."""
    return sum(g for g, _ in itertools.groupby(gzip.open(filename, 'rt'), key=is_header))


def count_fasta_files(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        return sum(g for g, _ in itertools.groupby(f, key=is_header))


def str_punctuation_strip(word):
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for _ in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


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

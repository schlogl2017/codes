#!usr/bin/env python
#-*- encoding: utf-8 -*-


import gzip
import itertools


def is_header(line):
    return line[0] == '>'

def parse_multi_fasta_file(filename):
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
        fasta_iter = (it[1] for it in itertools.groupby(f, is_header))
        for name in fasta_iter:
            name = next(name)
            sequences = ''.join(seq.strip() for seq in next(fasta_iter))
            yield name.strip()[1:], sequences


def str_punctuation_strip(word):
    """It strip all the punctuation and other characters in a string."""
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for char in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


# filename = 'genomes/refseq/test/GCF_000005825.2/GCF_000005825.2_ASM582v2_cds_from_genomic.fna.gz'
# filename = 'genomes/refseq/test/fastaFile.fa'
# for name, seq in parse_multi_fasta_file(filename):
#     name = '/'.join(str_punctuation_strip(name)[1:4:2])
#     print(name)
#     print(seq)


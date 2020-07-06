#!/usr/bin/env python
# coding: utf-8


import os
import gzip
import argparse
import itertools
from collections import defaultdict, Counter
import codons
from codons import codon_table_std_RNA


parser = argparse.ArgumentParser(description='Counting codons in genomes')
group = parser.add_mutually_exclusive_group()
parser.add_argument('--path', type=str, help='Path for the fasta files')
parser.add_argument('--output', type=str, help='Filename for the csv output file')
parser.add_argument('--length', type=int, help='Integer representing the mers length')
#parser.add_argument('--codon_map', type=str, help='Filename for the csv output file')
group.add_argument('--count', action='store_true', help='If opted the script will count the codons.')
group.add_argument('--list', action='store_true', help='If opted the script will only returns a list of codons.')
group.add_argument('--translation', action='store_true', help='If opted the script will only returns a list of proteins.')
args = parser.parse_args()


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


def str_punctuation_strip(word):
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for char in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


def get_translation(sequence, codon_map):
    lst_codons = codons.get_codons(sequence)
    protein = [codon_map[codon] for codon in lst_codons if codon not in ['UAA', 'UAG', 'UGA']]
    return ''.join(protein)


def get_kmers(sequence, k=1):
    """Yields all the kmers of length k from the sequence input.
    sequence = a string representing a DNA, Protein or RNA
    k = integer representing the length of the kmer/substrings    
    """
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        yield kmer


dir_name = args.path

infiles = []
for path, subdirs, files in os.walk(dir_name):
    for name in files:
        input_files = os.path.join(path, name)
        if input_files.endswith('.gz'):
            infiles.append(input_files)


cnt_files = 0
proteins = defaultdict(list)
kmers = defaultdict(list)
kmers_count = defaultdict(Counter)

for filename in infiles:
    print('Starting reading the fasta files', filename)
    for name, sequence in parse_multi_fasta_file_compressed_or_not(filename):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        proteins = get_translation(sequence, codon_table_std_RNA)
        kmers = get_kmers(proteins, args.length)
        if args.count:
            kmers_count[name].update(kmers)
        elif args.list:
            kmers.setdefault(name, [])
            kmers[name].append(', '.join(kmers(sequence, args.length)))
        elif args.translation:
            proteins[name].append(''.join(proteins))
    cnt_files += 1




#!usr/bin/env python

import os
import argparse
from collections import defaultdict, Counter
from toolz import partition, sliding_window
from toolz import itertoolz
import read_fasta_uniprot
import codons
from codons import codon_table_std_DNA, codon_table_std_RNA


def my_codons(sequence, mol_type='RNA'):
    """Return a generator of all codons(substring of length 3) with no-overlap window
    from a sequence(string) of DNA/RNA."""
    seq = sequence.upper()
    if mol_type == 'RNA':
        seq = seq.replace('T', 'U')
        return (''.join(c) for c in partition(3, seq))
    elif mol_type == 'DNA':
        return (''.join(c) for c in partition(3, seq))


def kmers(sequence, k):
    """Returns a generator of all mers(substring) of length k with overlap window
    from a string."""
    return (''.join(c) for c in sliding_window(k, sequence))


def count_subsatrings(sequence, k):
    """Returns the count of substrings of k length from a string"""
    mers = kmers(sequence, k)
    return itertoolz.frequencies(mers)


def get_translation(sequence, mol_type):
    """Translate the codons obtained from a fasta file with genomes CDS."""
    stop_codons = ['UAA', 'UAG', 'UGA', 'TAA', 'TAG', 'TGA']
    codon_list = my_codons(sequence, mol_type)
    if mol_type == 'RNA':
        codon_map = codon_table_std_RNA
    else:
        codon_map = codon_table_std_DNA
    protein = []
    for codon in codon_list:
        if codon not in codon_map:
            protein.append('?')
        elif codon in codon_map and codon not in stop_codons:
            protein.append(codon_map[codon])
    return ''.join(protein)


parser = argparse.ArgumentParser(description='Analysys of CDS genomes')
parser.add_argument('--path', type=str, help='Path for the directory of fasta files')
parser.add_argument('--output', type=str, help='Filename for the csv output file')
parser.add_argument('--length', type=int, help='Integer representing the mers/substring length')
parser.add_argument('--codon_table', type=str, help='Option of the codon table(DNA/RNA.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--count', action='store_true', help='Option for count the mers/substrings.')
group.add_argument('--list', action='store_true', help='Option to obtain the list of mers/substrings.')
group.add_argument('--translation', action='store_true', help='Option to obtain the list of '
                                                              'translated cds from a genome.')
args = parser.parse_args()


# find the path to the files
dir_name = args.path


# making a list of the fasta files
infiles = []
for path, subdirs, files in os.walk(dir_name):
    for name in files:
        input_files = os.path.join(path, name)
        if input_files.endswith('.gz'):
            infiles.append(input_files)

# put the sequences or mers /mers_counts in a dict
mers_count = defaultdict(Counter)
mers_lst = defaultdict(list)
prot_lst = defaultdict(list)
cnt_files = 0

for filename in infiles:
    print('Starting reading the fasta files', filename)
    for name, sequence in read_fasta_uniprot.parse_multi_fasta_file_compressed_or_not(filename):
        name = '/'.join(read_fasta_uniprot.str_punctuation_strip(name)[1:4:2])
        protein = get_translation(sequence, codons.codon_table_std_DNA)
        if args.count:
            mers_count[name].update(kmers(protein, args.length))
        elif args.list:
            mers_lst[name].append(list(kmers(protein, args.length)))
        elif args.translation:
            prot_lst[name].append(get_translation(sequence, args))
    cnt_files += 1


print('Writing the csv file.')
with open(args.path + '/' + args.output, 'w') as fout:
    if args.count:
        for key, count in mers_count.items():
            fout.write(key + ',' + str(count) + '\n')
    elif args.list:
        for key, lst in mers_lst.items():
            fout.write(key + ',' + str(lst) + '\n')
    elif args.translation:
        for key, lst in prot_lst.items():
            fout.write(key + ',' + str(lst) + '\n')
print('Were read and manipulated: {} fasta files'.format(cnt_files))
print('Job done')

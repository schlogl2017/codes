#!usr/bin/env python


import os
import argparse
from collections import defaultdict, Counter
from more_itertools import windowed
from toolz import partition
import read_fasta_uniprot
from codons import codon_table_std_RNA, codon_table_std_DNA


parser = argparse.ArgumentParser(description='Translate codons and get the mers from genomes')
parser.add_argument('--path', '-pat', type=str, required=True,
                    help='Path for the fasta files')
parser.add_argument('--output', '-out', type=str, required=True,
                    help='Filename for the csv output file')
parser.add_argument('--mol_type', type=str, required=True,
                    help='Path for the fasta files')
parser.add_argument('--length', '-len', type=int, dest='length',
                    help='Integer representing the mers length')
parser.add_argument('--fasta', '-faa', type=str, help='Return a fasta file')
group = parser.add_mutually_exclusive_group()
group.add_argument('--count', '-cnt', action='store_true',
                   help='If opted the script will count the mers.')
group.add_argument('--list', '-lt', action='store_true',
                   help='If opted the script will only returns a list of mer.')
group.add_argument('--translation', '-trans', action='store_true',
                   help='Return the list of translated cds from a genome.')
args = parser.parse_args()


def my_codons(sequence, mol_type='RNA'):
    """Return a generator of all codons(substring of length 3) with no-overlap window
    from a sequence(string) of DNA/RNA."""
    seq = sequence.upper()
    if mol_type == 'RNA':
        seq = seq.replace('T', 'U')
        return (''.join(c) for c in partition(3, seq))
    elif mol_type == 'DNA':
        return (''.join(c) for c in partition(3, seq))


# get the proteins from a fasta file with genomes cds
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


# get the mers in the fasta file with genomes cds
def get_mers(sequence, k=1):
    """returns the mers of length k obtained from a fasta file with genomes CDS or translated proteins."""
    return (''.join(mers) for mers in windowed(sequence, k))


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
protein_lst = defaultdict(list)
protein_dict = defaultdict()
cnt_files = 0

for filename in infiles:
    print('Starting reading the fasta file', filename)
    for name, sequence in read_fasta_uniprot.parse_multi_fasta_file_compressed_or_not(filename):
        name = '/'.join(read_fasta_uniprot.str_punctuation_strip(name)[1:4:2])
        protein = get_translation(sequence, args.mol_type)
        kmers = [mer for mer in get_mers(sequence, 5)]
        if args.count:
            mers_count[name].update(kmers)
        elif args.list:
            mers_lst[name].append(kmers)
        elif args.translation:
            protein_lst[name].append(protein)
        elif args.fasta:
            protein_dict[name].update(protein)
    cnt_files += 1

print('Writing the csv file.')
with open(args.path + '/' + args.output, 'w') as fout:
    if args.count:
        for key, val in mers_count.items():
            fout.write(key + ',' + str(val) + '\n')
# with open(args.path + '/' + args.output, 'w') as fout:
#     if args.count:
#         for key, count in mers_count.items():
#             fout.write(key + ',' + str(count) + '\n')
#     elif args.list:
#         for key, lst in mers_lst.items():
#             fout.write(key + ',' + str(lst) + '\n')
#     elif args.translation:
#         for key, lst in protein_lst.items():
#             fout.write(key + ',' + str(lst) + '\n')
#     elif args.fasta:
#         for key, seq in protein_lst.items():
#             fout.write('>' + key + '\n' + str(seq) + '\n')


print('Were read and manipulated: {} fasta files'.format(cnt_files))
print('Job done')

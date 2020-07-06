#!usr/bin/env python


import os
import argparse
from collections import defaultdict, Counter
from toolz import partition
from toolz import itertoolz
import read_fasta_uniprot
from codons import codon_table_std_RNA, codon_table_std_DNA


def my_codons(sequence, mol_type='RNA'):
    """Return a generator of all codons(substring of length 3) with no-overlap window
    from a sequence(string) of DNA/RNA."""
    seq = sequence.upper()
    if mol_type == 'RNA':
        seq = seq.replace('T', 'U')
        return (''.join(c) for c in partition(3, seq))
    elif mol_type == 'DNA':
        return (''.join(c) for c in partition(3, seq))


def count_codons(sequence, mol_type):
    """Returns the count of codons from a fasta file with CDS from a genome."""
    codons = my_codons(sequence, mol_type)
    return itertoolz.frequencies(codons)


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
parser.add_argument('--path', '-p', type=str, required=True,
                    help='Path for the directory of fasta files')
parser.add_argument('--output', '-out', type=str, required=True,
                    help='Filename for the csv output file')
parser.add_argument('--mol_type', type=str, dest='mol_type', required=True,
                    help='Option of the codon table(DNA/RNA.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--count', '-cnt', action='store_true',
                   help='Option for count the mers/substrings.')
group.add_argument('--list', '-lst', action='store_true', dest='list',
                   help='Option to obtain the list of mers/substrings.')
group.add_argument('--translation', '-trans', action='store_true',
                   help='Option to obtain the list of translated cds from a genome.')
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
        elif input_files.endswith('.fna') or input_files.endswith('.faa'):
            infiles.append(input_files)
        else:
            nozip_infiles.append(input_files)


# put the sequences or mers /mers_counts in a dict
codons_count = defaultdict(Counter)
codons_lst = defaultdict(list)
prot_lst = defaultdict(list)
cnt_files = 0

for filename in infiles:
    print('Starting reading the fasta file: ', filename)
    for name, sequence in read_fasta_uniprot.parse_multi_fasta_file_compressed_or_not(filename):
        name = '/'.join(read_fasta_uniprot.str_punctuation_strip(name)[1:4:2])
        protein = get_translation(sequence, args.mol_type)
        codons = my_codons(sequence, args.mol_type)
        if args.count:
            codons_count[name].update(codons)
        elif args.list:
            codons_lst[name].append(codons)
        elif args.translation:
            prot_lst[name].append(protein)
    cnt_files += 1


print('Writing the csv file.')
with open(args.path + '/' + args.output, 'w') as fout:
    if args.count:
        for key, count in codons_count.items():
            fout.write(key + ',' + str(count) + '\n')
    elif args.list:
        for key, lst in codons_lst.items():
            fout.write(key + ',' + str(lst) + '\n')
    elif args.translation:
        for key, lst in prot_lst.items():
            fout.write(key + ',' + str(lst) + '\n')


print('Were read and manipulated: {} fasta files.'.format(cnt_files))
print('Job done.')

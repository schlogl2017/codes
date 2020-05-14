#!usr/bin/env python


import os
import argparse
import time
from toolz import partition
from collections import defaultdict, Counter
from fasta_reader import parse_multi_fasta_file, str_punctuation_strip
from std_codons_tables import codon_table_std_RNA, codon_table_std_DNA, stop_codons_std


start_time = time.process_time()


parser = argparse.ArgumentParser(description='Kmers from proteins.')
parser.add_argument('--path', '-p', type=str, required=True, dest='path', help='Path to the files')
parser.add_argument('--out', '-o', type=str, required=True, dest='output', help='name output files.')
parser.add_argument('--moltype', '-mt', type=str, required=True, dest='mol_type', help='Type mol DNA/RNA')
parser.add_argument('--length', '-l', type=int, dest='k', help='Length of the kmers')
group = parser.add_mutually_exclusive_group()
group.add_argument('--count', '-c', action='store_true', help='Count kmers')
group.add_argument('--list', '-ls', action='store_true', help='List of kmer')
group.add_argument('--fasta', '-faa', action='store_true', help='If you need a protein fasta file.')
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
def get_translation(sequence, mol_type='RNA'):
    """Translate the codons obtained from a fasta file with genomes CDS."""
    codon_list = my_codons(sequence, mol_type)
    if mol_type == 'RNA':
        codon_map = codon_table_std_RNA
    else:
        codon_map = codon_table_std_DNA
    protein = []
    for codon in codon_list:
        if codon not in codon_map:
            protein.append('?')
        elif codon in codon_map and codon not in stop_codons_std:
            protein.append(codon_map[codon])
    return ''.join(protein)


def get_kmers(sequence, k):
    """Yields all the kmers of length k from the sequence input.
    sequence = a string representing a DNA, Protein or RNA
    k = integer representing the length of the kmer/substrings
    """
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        yield kmer


def get_list_files(dirName):
    """Retruns a list of file and sub directories names in thget_list_filese given directory"""
    list_files = os.listdir(dirName)
    all_files = list()
    # Iterate over all the entries
    for entry in list_files:
        # Create full path
        full_path = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files = all_files + get_list_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


# making a generator of the fasta files
def check_filenames(files):
    zipped = []
    not_zip = []
    for file in files:
        if file.endswith('.gz'):
            zipped.append(file)
        elif file.endswith('.fa') or file.endswith('.faa'):
            not_zip.append(file)
    return zipped, not_zip


dir_name = args.path
k = args.k

# getting all the files
all_files = get_list_files(dir_name)

print('Start reading the fasta files from directories: {}'.format(dir_name))
files_z, _ = check_filenames(all_files)

prot_mers_counter = Counter()
prot_mers_lst = defaultdict(list)
prot_dict = defaultdict()
cnt = 0
for file in files_z:
    print('Fasta file: {}'.format(file))
    for name, sequence in parse_multi_fasta_file(file):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        proteins = get_translation(sequence, args.mol_type)
        if args.count:
            prot_mers_counter.update(list(mer for mer in get_kmers(proteins, args.k)))
        elif args.list:
            prot_mers_lst[name] = [mer for mer in get_kmers(proteins, args.k)]
        elif args.fasta:
            prot_dict[name] = proteins
    cnt += 1


print('Start writing the {} file.'.format(args.output))
with open(dir_name + '/' + args.output, 'w') as fout:
    if args.count:
        for (mers, count) in prot_mers_counter.items():
            fout.write(mers + ',' + str(count) + '\n')
    elif args.list:
        for (name, mers) in prot_mers_lst.items():
            fout.write(name + ',' + str(mers) + '\n')
    elif args.fasta:
        for (name, prot) in prot_dict.items():
            fout.write('>' + name + '\n' + prot + '\n')


print('Were read and manipulated: {} fasta files'.format(cnt))
print('Time of execution: {} seconds'.format(time.process_time() - start_time))
print('Script finished!')

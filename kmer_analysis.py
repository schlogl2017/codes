#!usr/bin/env python


import os
import argparse
from collections import Counter, defaultdict
from toolz import sliding_window
from read_fasta_uniprot import parse_multi_fasta_file_compressed_or_not
from read_fasta_uniprot import str_punctuation_strip

parser = argparse.ArgumentParser(description='Count the mers in CDS genomes')
parser.add_argument('-k', '--length', type=int, help='Integer representing th mers length.')
parser.add_argument('-o', '--output', type=str, help='Output filename.')
parser.add_argument('-pa', '--path', type=str, help='Path for the directory with the files.')
group = parser.add_mutually_exclusive_group()
group.add_argument('-c', '--count', action='store_true', help='Count the mers of k length.')
group.add_argument('-lst', '--list', action='store_true', help='Return a list of mers of k length')
args = parser.parse_args()


def kmers(sequence, k):
    """Returns a generator of all mers(substring) of length k with overlap window
    from a string."""
    return (''.join(c) for c in sliding_window(k, sequence))


# path for the directory with the fastas files
dir_name = args.path


# making a list of the fasta files
infiles = []
for path, subdirs, files in os.walk(dir_name):
    for name in files:
        input_files = os.path.join(path, name)
        if input_files.endswith('.gz'):
            infiles.append(input_files)


# analysis
count = defaultdict(Counter)
num_files = 0
for file in infiles:
    print('Reading filename: {}'.format(file))
    for name, seq in parse_multi_fasta_file_compressed_or_not(file):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        count[name].update(kmers(seq, args.length))
    num_files += 1


print('Were processed {} files'.format(num_files))

with open(args.path + '/' + args.output, 'w') as fhout:
    for k, v in count.items():
        fhout.write(k + ',' + str(v) + '\n')

# checking if the file was created
filename = os.path.isfile(args.path + '/' + args.output + '.csv')
print(filename)

print('Job finished.')

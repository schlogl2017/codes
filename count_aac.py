#!usr/bin/env python

import gzip
from collections import defaultdict, Counter
from Bio import SeqIO





def count_aac_id_sequence(filename):
    """Yield a dict of dict with the count of aminoacids from a compressed fasta file. 
    The keys are the name of the sequences and the values as the aac counts."""
    with gzip.open(filename, "rt") as handle:
        aac_count = defaultdict(Counter)
        for record in SeqIO.parse(handle, "fasta"):
            aac_count[record.id].update(record.seq)
            yield aac_count

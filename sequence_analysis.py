#!usr/bin/env python

from toolz import partition, sliding_window
from toolz import itertoolz
from collections import Counter


def codons(sequence, mol_type='RNA'):
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


def count_stuffs(sequence, k):
    mers = kmers(sequence, k)
    return itertoolz.frequencies(mers)









codon_table_std_RNA = {'AAA': 'K',
                       'AAC': 'N',
                       'AAG': 'K',
                       'AAU': 'N',
                       'ACA': 'T',
                       'ACC': 'T',
                       'ACG': 'T',
                       'ACU': 'T',
                       'AGA': 'R',
                       'AGC': 'S',
                       'AGG': 'R',
                       'AGU': 'S',
                       'AUA': 'I',
                       'AUC': 'I',
                       'AUG': 'M',
                       'AUU': 'I',
                       'CAA': 'Q',
                       'CAC': 'H',
                       'CAG': 'Q',
                       'CAU': 'H',
                       'CCA': 'P',
                       'CCC': 'P',
                       'CCG': 'P',
                       'CCU': 'P',
                       'CGA': 'R',
                       'CGC': 'R',
                       'CGG': 'R',
                       'CGU': 'R',
                       'CUA': 'L',
                       'CUC': 'L',
                       'CUG': 'L',
                       'CUU': 'L',
                       'GAA': 'E',
                       'GAC': 'D',
                       'GAG': 'E',
                       'GAU': 'D',
                       'GCA': 'A',
                       'GCC': 'A',
                       'GCG': 'A',
                       'GCU': 'A',
                       'GGA': 'G',
                       'GGC': 'G',
                       'GGG': 'G',
                       'GGU': 'G',
                       'GUA': 'V',
                       'GUC': 'V',
                       'GUG': 'V',
                       'GUU': 'V',
                       'UAA': '*',
                       'UAC': 'Y',
                       'UAG': '*',
                       'UAU': 'Y',
                       'UCA': 'S',
                       'UCC': 'S',
                       'UCG': 'S',
                       'UCU': 'S',
                       'UGA': '*',
                       'UGC': 'C',
                       'UGG': 'W',
                       'UGU': 'C',
                       'UUA': 'L',
                       'UUC': 'F',
                       'UUG': 'L',
                       'UUU': 'F'}


codon_table_std_DNA = {'AAA': 'K',
                       'AAC': 'N',
                       'AAG': 'K',
                       'AAT': 'N',
                       'ACA': 'T',
                       'ACC': 'T',
                       'ACG': 'T',
                       'ACT': 'T',
                       'AGA': 'R',
                       'AGC': 'S',
                       'AGG': 'R',
                       'AGT': 'S',
                       'ATA': 'I',
                       'ATC': 'I',
                       'ATG': 'M',
                       'ATT': 'I',
                       'CAA': 'Q',
                       'CAC': 'H',
                       'CAG': 'Q',
                       'CAT': 'H',
                       'CCA': 'P',
                       'CCC': 'P',
                       'CCG': 'P',
                       'CCT': 'P',
                       'CGA': 'R',
                       'CGC': 'R',
                       'CGG': 'R',
                       'CGT': 'R',
                       'CTA': 'L',
                       'CTC': 'L',
                       'CTG': 'L',
                       'CTT': 'L',
                       'GAA': 'E',
                       'GAC': 'D',
                       'GAG': 'E',
                       'GAT': 'D',
                       'GCA': 'A',
                       'GCC': 'A',
                       'GCG': 'A',
                       'GCT': 'A',
                       'GGA': 'G',
                       'GGC': 'G',
                       'GGG': 'G',
                       'GGT': 'G',
                       'GTA': 'V',
                       'GTC': 'V',
                       'GTG': 'V',
                       'GTT': 'V',
                       'TAA': '*',
                       'TAC': 'Y',
                       'TAG': '*',
                       'TAT': 'Y',
                       'TCA': 'S',
                       'TCC': 'S',
                       'TCG': 'S',
                       'TCT': 'S',
                       'TGA': '*',
                       'TGC': 'C',
                       'TGG': 'W',
                       'TGT': 'C',
                       'TTA': 'L',
                       'TTC': 'F',
                       'TTG': 'L',
                       'TTT': 'F'}


def get_all_codons(sequence, mol_type='RNA'):
    """ Returns a list of triplets that represent all the possible codons in a DNA or RNA sequence. """

    seq = sequence.upper()
    codons = []
    if mol_type == 'RNA':
        seq = seq.replace('T', 'U')
    elif mol_type == 'DNA':
        seq = seq

    # sliding window of 1
    for i, char in enumerate(seq):
        codon = seq[i:i+3]
        if (mol_type == 'RNA') and (codon in codon_table_std_RNA) and (len(codon) == 3):
            codons.append(codon)
        else:
            if (mol_type == 'DNA') and (codon in codon_table_std_DNA) and (len(codon) == 3):
                codons.append(codon)
    return codons


def get_codons(sequence, mol_type='RNA'):
    """ Returns a list of triplets that represent all the possible codons in a DNA or RNA sequence. without
    a sliding window"""

    seq = sequence.upper()
    codons = []
    if mol_type == 'RNA':
        seq = seq.replace('T', 'U')
    elif mol_type == 'DNA':
        seq = seq

    for i in range(0, (len(seq) - 2), 3):
        codon = seq[i:i+3]
        if mol_type == 'RNA' and codon in codon_table_std_RNA:
            codons.append(codon)
        else:
            if mol_type == 'DNA' and codon in codon_table_std_DNA:
                codons.append(codon)
    return codons


def test():
    assert get_all_codons('ACGTAGCCAATGCGCCAA', 'RNA') == ['ACG', 'CGU', 'GUA',
                                                           'UAG', 'AGC', 'GCC',
                                                           'CCA', 'CAA', 'AAU',
                                                           'AUG', 'UGC', 'GCG',
                                                           'CGC', 'GCC', 'CCA',
                                                           'CAA'], "test failed"
    assert get_all_codons('ACGTAGCCAATGCGCCAA', 'DNA') == ['ACG', 'CGT', 'GTA',
                                                           'TAG', 'AGC', 'GCC',
                                                           'CCA', 'CAA', 'AAT',
                                                           'ATG', 'TGC', 'GCG',
                                                           'CGC', 'GCC', 'CCA',
                                                           'CAA'], "test failed"
    assert get_all_codons('acgtagccaatgcgccaa', 'RNA') == ['ACG', 'CGU', 'GUA',
                                                           'UAG', 'AGC', 'GCC',
                                                           'CCA', 'CAA', 'AAU',
                                                           'AUG', 'UGC', 'GCG',
                                                           'CGC', 'GCC', 'CCA',
                                                           'CAA'], "test failed"
    assert get_all_codons('acgtagccaatgcgccaa', 'DNA') == ['ACG', 'CGT', 'GTA',
                                                           'TAG', 'AGC', 'GCC',
                                                           'CCA', 'CAA', 'AAT',
                                                           'ATG', 'TGC', 'GCG',
                                                           'CGC', 'GCC', 'CCA',
                                                           'CAA'], "test failed"
    assert get_all_codons('', 'DNA') == [], "test failed"
    assert get_all_codons('', 'RNA') == [], "test failed"
    assert get_codons('', 'DNA') == [], "test failed"
    assert get_codons('', 'RNA') == [], "test failed"
    assert get_all_codons('ACGTGANN', 'DNA') == ['ACG', 'CGT', 'GTG', 'TGA'], "test failed"
    assert get_all_codons('AUGVVAA', 'RNA') == ['AUG'], "test failed"
    assert get_all_codons('VVAA', 'RNA') == [], "test failed"
    assert get_codons('ACGTGANN', 'DNA') == ['ACG', 'TGA'], "test failed"
    assert get_codons('AUGVVAA', 'RNA') == ['AUG'], "test failed"
    assert get_codons('VAA', 'RNA') == [], "test failed"
    assert get_codons('VAA', 'DNA') == [], "test failed"
    assert get_all_codons('VAA', 'RNA') == [], "test failed"
    assert get_all_codons('VAA', 'DNA') == [], "test failed"

def main():
    test()


if __name__ == '__main__':
    main()
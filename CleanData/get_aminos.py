import os
import sys
sys.path.append('../CUB_Code')

import findSequenceById as FSBID
from CodonLibraries import loadSequence

codonTable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

seq_path = os.path.join('..', 'Data', 'fnaCollection')
num_path = os.path.join('..', 'Data', 'fna_txt')

seq_list = [os.path.join(seq_path, f) for f in os.listdir(seq_path)]

all_seqs = []
for f in seq_list:
    seq_tups = [loadSequence(l[1]) for l in FSBID.findSequenceByID(f)]
    all_seqs += seq_tups

print(len(all_seqs))

amino_seqs = []
for s in all_seqs:
    amino_seqs.append([codonTable[si] for si in s])

print(len(amino_seqs))
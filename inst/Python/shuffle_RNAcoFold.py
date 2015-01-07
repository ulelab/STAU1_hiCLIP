#!/usr/bin/python 
# This script is downloaded from http://code.google.com/p/shortread-toolbox/source/browse/trunk/scripts/shuffle-fasta.py?r=8 and minor modifications are added.

import sys
import os
import random

random.seed(1)

usage_mesg = 'Usage: shuffle-fasta.py <fasta file>'
if( len(sys.argv) != 3 ):
    print usage_mesg
    sys.exit(1)

filename_fasta = sys.argv[1]
if( not os.access(filename_fasta,os.R_OK) ):
    print "%s is not accessible"%filename_fasta
    print usage_mesg
    sys.exit(1)

filename_shuffle = sys.argv[2]

f_fasta = open(filename_fasta,'r')
header = ''
seq = dict()

for line in f_fasta:
    if(line.startswith('>')):
        tokens = line.strip().split()
        # print tokens
        id_tokens = tokens[0].lstrip('>')
        # print id_tokens
        header = ''+id_tokens
        if header in seq.keys():
            raise ValueError("Duplicated seq names")
        seq[header] = ''
    else:
        seq[header] += line.strip()

f_fasta.close()

def shuffle_seq(rna_seq):
    seq_list = list(rna_seq)
    random.shuffle(seq_list)
    random_rnaseq = ''.join(seq_list)
    return random_rnaseq

f_shuffle = open(filename_shuffle,'w')

for h in seq.keys():
    seqs = seq[h].split('&')
    left_random = shuffle_seq(seqs[0])
    right_random = shuffle_seq(seqs[1])
    random_seq = '&'.join((left_random, right_random))

    f_shuffle.write(">%s\n%s\n"%(h,random_seq))

# for h in seq.keys():
#     seqs = seq[h]
#     left_random = shuffle_seq(seqs)
#     random_seq = left_random
#
#     f_shuffle.write(">%s\n%s\n"%(h,random_seq))

f_shuffle.close()
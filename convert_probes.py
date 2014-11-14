#!/usr/bin/env python

'''
Convert Illumina probes to fastq format.

Example usage:
 python convert_probes.py data/HumanHT-12_V4_0_R2_15002873_B.txt > data/ht12_probes.fq

Example output:

@NM_017583.3:TRIM44
CCTGCCTGTCTGCCTGTGACCTGTGTACGTATTACAGGCTTTAGGACCAG
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

import sys

file = sys.argv[1]
manifest = open(file, 'r')

for line in manifest:
    cols = line.strip().split('\t')
    if cols[0] == 'Homo sapiens':
        name = '@' + cols[13] + ':' + cols[2] + ':' + cols[11]
        seq = cols[17]
        qual = '~' * len(seq)
        sys.stdout.write(name + '\n' + seq + '\n' + '+\n' + qual + '\n')

manifest.close()

#!/usr/bin/env python

'''
Convert HapMap SNPs to bed format.
Include reference allele frequency as score.

Example usage:
  python convert_snps.py data/genotypes_chrX_CEU_r28_nr.b36_fwd.txt.gz  data/allele_freqs_chrX_CEU_r28_nr.b36_fwd.txt.gz > data/chrX_snps.bed

Example output:
chrX	90117	90118	rs6608427	+	0.415
chrX	94674	94675	rs6608419	+	0.423
chrX	94752	94753	rs6608418	+	0.413
chrX	96873	96874	rs6608411	+	1.000
chrX	108566	108567	rs6655866	+	1.000
chrX	108821	108822	rs6649842	+	1.000
chrX	109804	109805	rs6423165	+	0.420
chrX	109921	109922	rs6608381	+	0.411
'''

import sys
import gzip

# Open connection and skip header line
geno = gzip.open(sys.argv[1], 'r')
geno.readline()
freq = gzip.open(sys.argv[2], 'r')
freq.readline()

# Read in allele frequencies
d_freq = {}
for f in freq:
    f = str(f, encoding='utf8')
    f_cols = f.strip().split(' ')
    d_freq[f_cols[0]] = f_cols[11]

# Read in genotypes and write out in bed format with allele frequency
for g in geno:
    g = str(g, encoding='utf8')
    g_cols = g.strip().split(' ')
    chrom = g_cols[2]
    start = str(int(g_cols[3]) - 1)
    end = g_cols[3]
    name = g_cols[0]
    strand = '+'
    score = d_freq.get(name, 'NA')
    sys.stdout.write(chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + score + '\t' + strand + '\n')

geno.close()
freq.close()

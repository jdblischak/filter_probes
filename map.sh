#!/bin/bash

# Map probes to genome.

# Index genome.
software/bwa/bwa index data/hg19.fa

# Convert probes to fastq format.
python convert_probes.py data/HumanHT-12_V4_0_R2_15002873_B.txt > data/ht12_probes.fq

# Map with bwa (use backtrack algorithm because reads are 50 bp)
software/bwa/bwa aln data/hg19.fa data/ht12_probes.fq > data/ht12_probes.sai
software/bwa/bwa samse data/hg19.fa data/ht12_probes.sai data/ht12_probes.fq > data/ht12_probes.sam

# Filter with samtools
samtools view -SbF 4 data/ht12_probes.sam > data/ht12_probes.bam

# Convert to bed format
bedtools bamtobed -i data/ht12_probes.bam > data/ht12_probes.bed

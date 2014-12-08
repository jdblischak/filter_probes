## Purpose

Single nucleotide polymorphisms (SNPs) common in a population can
affect the hybridation kinetics of gene expression microarrays. This
is especially troublesome for studies that investigate differences in
gene expression between individuals.

## Technical description

For our analysis, we used a subset of the 47,321 probes designed to
target genome of the Illumina HT-12v4 bead array that we determined
were not affected by common SNPs in the HapMap CEU population. First,
we mapped the probes to human genome hg19 and kept only those reads
which matched perfectly (39,801 probes; note that this also resulted
in the removal of any probes spanning exon-exon junctions). Second, we
downloaded the HapMap CEU SNPs
(http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2010-08_phaseII+III/forward/)
and converted their coordinates from hg18 to hg19 using the UCSC
liftOver utility. We kept only those probes which did not overlap any
SNP with a minor allele frequency of at least 5% (36,893
probes). Third, we converted the Illumina probe IDs to Ensembl gene
IDs using the R/Bioconductor package biomaRt ([Durinck et al.,
2009][Durinck2009]) and kept only those probes which are associated
with exactly one Ensembl gene ID (Ensembl 75 - Feb 2014; 23,890 probes). The full pipeline
was implemented using the Python package Snakemake ([Koster & Rahmann,
2012][Koster2012]).

[Durinck2009]: http://www.nature.com/nprot/journal/v4/n8/full/nprot.2009.97.html
[Koster2012]: http://bioinformatics.oxfordjournals.org/content/28/19/2520.long

## Results

* Number of probes designed for human genome: 47231 
* Number of probes mapped to hg19: 43143 
* Number of probes mapped uniquely to hg19: 39801 
* Number of probes mapped with quality score >= 37 to hg19: 37963 
* Number of probes mapped without a SNP ( MAF >= 0.05 ): 42027 
* Number of probes mapped with quality score >= 37 to hg19 and without a SNP ( MAF >= 0.05 ): 36893 
* Number of probes associated with a unique Ensembl gene ID: 23890

## Software

This analysis pipeline requires the following to be installed:

* Python3
* Snakemake
* Samtools
* bedtools

The pipeline installs the following:

* bwa
* liftOver

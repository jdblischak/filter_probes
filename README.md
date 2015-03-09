## Purpose

Single nucleotide polymorphisms (SNPs) common in a population can
affect the hybridation kinetics of gene expression microarrays. This
is especially troublesome for studies that investigate differences in
gene expression between individuals.

## Usage

Clone the repository:

```
git clone git@github.com:jdblischak/filter_probes.git
cd filter_probes
```

To run sequentially, run:

```
snakemake -s make.py
```

To run in parallel, use the `-c` and `-j` options to submit each job
to your grid computing system. This will vary based on your
setup. Here is an example using grid engine on my work cluster:

```
snakemake -s make.py -j 30 -c "qsub -l h_vmem={params.h_vmem} -N {params.name} -V -j y -cwd -o {log}"
```

This will submit up to 30 jobs at once, each requesting the virtual
memory and name specified for each job in the Snakefile. Since this
can take hours to complete, I recommend submitting the job in the
background immune to hangups, e.g.:

```
nohup snakemake -s make.py -j 30 -c "qsub -l h_vmem={params.h_vmem} -N {params.name} -V -j y -cwd -o {log}" &
```

## Technical description

For our analysis, we used a subset of the 47,321 probes designed to
target genome of the Illumina HT-12v4 bead array that we determined
were not affected by common SNPs in the European population. First,
we mapped the probes to human genome hg19 and kept only those with a
quality score of 37 (40,198 probes; note that we also explicitly
pre-filtered the 5,587 probes which were annotated as spanning
exon-exon junctions to avoid mapping errors). Second, we downloaded
the phase 1 integrated call sets from the 1000 Genomes Project
(ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/).
We kept only those probes which did not overlap any
SNP with a minor allele frequency in the European population of at least 5% (33,494
probes). Third, we converted the Illumina probe IDs to Ensembl gene
IDs using the R/Bioconductor package biomaRt ([Durinck et al.,
2009][Durinck2009]) and kept only those probes which are associated
with exactly one Ensembl gene ID (Ensembl 75 - Feb 2014; 21,676
probes). The full pipeline was implemented using the Python package
Snakemake ([Koster & Rahmann, 2012][Koster2012]).

[Durinck2009]: http://www.nature.com/nprot/journal/v4/n8/full/nprot.2009.97.html
[Koster2012]: http://bioinformatics.oxfordjournals.org/content/28/19/2520.long

## Results

* Number of probes designed for human genome: 47231
* Number of probes mapped to hg19: 40198
* Number of probes mapped uniquely to hg19: 37157
* Number of probes mapped with quality score >= 37 to hg19: 35529
* Number of probes mapped without a SNP ( MAF >= 0.05 ): 37903
* Number of probes mapped with quality score >= 37 to hg19 and without a SNP ( MAF >= 0.05 ): 33494
* Number of probes associated with a unique Ensembl gene ID: 21676

## Software

This analysis pipeline requires the following to be installed:

* Python3
* Snakemake
* Samtools
* bedtools
* R package argparse
* R/Bioconductor package biomaRt
* vcftools

The pipeline installs the following:

* bwa


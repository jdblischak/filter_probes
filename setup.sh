#!/bin/bash

# Install:
#   -bwa
#   -liftOver

# Download:
#   -HapMap CEU genotypes
#   -HapMap CEU allele frequencies
#   -Human genome hg19
#   -Illumina HT12 array manifest

# Installations
mkdir software
cd software

# Install bwa
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

# Install liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

cd ..

# Download data
mkdir data
cd data

# Download HapMap CEU genotypes
wget http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2010-08_phaseII+III/forward/genotypes_chr{1..22}_CEU_r28_nr.b36_fwd.txt.gz
wget http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2010-08_phaseII+III/forward/genotypes_chr{X,Y,M}_CEU_r28_nr.b36_fwd.txt.gz

# Download HapMap CEU allele frequencies
wget http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III/allele_freqs_chr{1..22}_CEU_r28_nr.b36_fwd.txt.gz
wget http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III/allele_freqs_chr{X,Y,M}_CEU_r28_nr.b36_fwd.txt.gz

# Download human genome hg19
rsync -avzuP globus.opensciencedatacloud.org::public/illumina/igenomes/Homo_sapiens/UCSC/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa .
mv genome.fa hg19.fa

# Download Illumina array manifest
wget http://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/humanht-12_v4_0_r2_15002873_b.txt.zip
unzip humanht-12_v4_0_r2_15002873_b.txt.zip

cd ..

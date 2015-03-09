

'''
To submit:
nohup snakemake -kps make.py -j 96 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log}" &
'''

import os
import shutil

# Paths
DATA_DIR = 'data/'
SRC_DIR = 'software/'
LOG_DIR = 'log/'

for d in [DATA_DIR, SRC_DIR, LOG_DIR]:
    if not os.path.isdir(d):
        os.mkdir(d)

# Settings
CHROM = [str(x) for x in range(1, 23)] + ['X', 'Y', 'M']
MAF = 0.05 # Minor Allele Frequency cutoff
MAP_SCORE = 37 # Mapping quality score cutoff

# Target rules
localrules: all

rule all:
	input: DATA_DIR + 'ht12_probes_snps_ceu_hg19_af_' + str(MAF) + '_map_' + str(MAP_SCORE) + '.txt'

rule test:
	input: DATA_DIR + 'genotypes_chr22.bed'

# Workflow
rule download_genos:
	output: DATA_DIR + 'ALL.chr{CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'download_genos' + wildcards.CHR
	log: LOG_DIR
	shell: 'wget -O {output} ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr{wildcards.CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz'

rule download_genos_Y_and_MT:
	output: DATA_DIR + 'ALL.chr{CHR}.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'download_genos' + wildcards.CHR
	log: LOG_DIR
	shell: 'wget -O {output} ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr{wildcards.CHR}.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz'

rule rename_files:
	input: expand(DATA_DIR + 'ALL.chr{CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz', CHR = CHROM[:23]),
               expand(DATA_DIR + 'ALL.chr{CHR}.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz', CHR = ['Y', 'MT'])
	output: expand(DATA_DIR + 'genotypes_chr{CHR}.vcf.gz', CHR = CHROM)
	params: h_vmem = '8g', bigio = '0',
                name = 'rename_files'
	log: LOG_DIR
	run: 
          for i in range(len(CHROM)):
             shutil.copy(input[i], output[i])

rule extract_geno_and_af:
	input: DATA_DIR + 'genotypes_chr{CHR}.vcf.gz'
	output: DATA_DIR + 'genotypes_chr{CHR}.txt.gz'
	params: h_vmem = '8g', bigio = '0',
                name = lambda wildcards: 'extract_geno_and_af.' + wildcards.CHR
	log: LOG_DIR
	shell: 'vcftools --gzvcf {input} --get-INFO AF --get-INFO EUR_AF \
                --get-INFO VT --stdout | gzip -c > {output}'

rule geno_to_bed:
	input: DATA_DIR + 'genotypes_chr{CHR}.txt.gz'
	output: DATA_DIR + 'genotypes_chr{CHR}.bed'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'geno_to_bed.' + wildcards.CHR
	log: LOG_DIR
	run:
          import gzip

          # Open connection and skip header line
          geno = gzip.open(input[0], 'r')
          geno.readline()

          # Open connection to output file
          bed = open(output[0], 'w')

          # VCF files for Y and M are completely different from
          # the other chromosomes. Writing separate functions
          # to process them.
          def process_chr(g_cols, bed):
              # Output as bed format.
              # g_cols: list of columns from input file
              # bed: handle to bed file

              eur_af = g_cols[5]
              variant_type = g_cols[6]
              if eur_af == '?' or variant_type != 'SNP':
                  return None
              chrom = 'chr' + g_cols[0]
              start = str(int(g_cols[1]) - 1)
              end = g_cols[1]
              name = chrom + '.' + end
              score = g_cols[5]
              strand = '+'
              bed.write(chrom + '\t' + start + '\t' + end + '\t' + \
                        name + '\t' + score + '\t' + strand + '\n')

          def process_chr_Y_MT(g_cols, bed):
              # Output as bed format.
              # Only for Y and MT b/c they do not have allele frequency
              # or variant type information.
              # g_cols: list of columns from input file
              # bed: handle to bed file

              # Discard if not SNP
              if len(g_cols[2]) != 1 or len(g_cols[3]) != 1:
                  return None
              # If MT, replace with M b/c that is name used by
              # hg19.
              if g_cols[0] == 'MT':
                  chrom = 'chrM'
              else:
                  chrom = 'chr' + g_cols[0]
              start = str(int(g_cols[1]) - 1)
              end = g_cols[1]
              name = chrom + '.' + end
              score = str(0.5) # Assign them all allele frequency of 0.5
                          # to be overly conservative.
              strand = '+'
              bed.write(chrom + '\t' + start + '\t' + end + '\t' + \
                        name + '\t' + score + '\t' + strand + '\n')
              
          # Read in genotypes and write out in bed format
          for g in geno:
              g = str(g, encoding='utf8')
              g_cols = g.strip().split('\t')
              if g_cols[0] in ['Y', 'MT']:
                  process_chr_Y_MT(g_cols, bed)
              else:
                  process_chr(g_cols, bed)

          geno.close()
          bed.close()

rule combine_genos:
	input: expand(DATA_DIR + 'genotypes_chr{CHR}.bed', CHR = CHROM)
	output: DATA_DIR + 'snps_EUR_1KG.bed'
	params: h_vmem = '8g', bigio = '0',
            name = 'combine_genos'
	log: LOG_DIR
	shell: 'cat {input} > {output}'

################################################################################
# Process probes
################################################################################

rule download_probes:
	output: DATA_DIR + 'humanht-12_v4_0_r2_15002873_b.txt.zip'
	params: h_vmem = '8g', bigio = '0',
                name = 'download_probes'
	log: LOG_DIR
	shell: 'wget -O {output} http://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/humanht-12_v4_0_r2_15002873_b.txt.zip'

rule unzip_probes:
	input: DATA_DIR + 'humanht-12_v4_0_r2_15002873_b.txt.zip'
	output: DATA_DIR + 'HumanHT-12_V4_0_R2_15002873_B.txt'
	params: h_vmem = '8g', bigio = '0',
                name = 'unzip_probes'
	log: LOG_DIR
	shell: 'unzip -p {input} > {output}'

'''
Convert Illumina probes to fastq format.

Example output:

@NM_017583.3:TRIM44
CCTGCCTGTCTGCCTGTGACCTGTGTACGTATTACAGGCTTTAGGACCAG
+
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''

rule convert_probes:
	input: DATA_DIR + 'HumanHT-12_V4_0_R2_15002873_B.txt'
	output: DATA_DIR + 'ht12_probes.fq'
	params: h_vmem = '8g', bigio = '0',
                name = 'convert_probes'
	log: LOG_DIR
	run:
          manifest = open(input[0], 'r')
          fastq = open(output[0], 'w')

          for line in manifest:
              cols = line.strip().split('\t')
              if cols[0] == 'Homo sapiens' and ':' not in cols[20]:
                  name = '@' + cols[13] + ':' + cols[2] + ':' + cols[8] + ':' + cols[11]
                  seq = cols[17]
                  qual = '~' * len(seq)
                  fastq.write(name + '\n' + seq + '\n' + '+\n' + qual + '\n')

          manifest.close()
          fastq.close()

rule install_bwa:
	output: SRC_DIR + 'bwa/bwa',
	params: h_vmem = '8g', bigio = '0',
            name = 'install_bwa'
	log: LOG_DIR
	shell: '''
        cd {SRC_DIR}
        git clone https://github.com/lh3/bwa.git
        cd bwa
        make
        '''

rule download_genome:
	output: DATA_DIR + 'hg19.fa',
	params: h_vmem = '8g', bigio = '0',
            name = 'download_genome'
	log: LOG_DIR
	shell: '''
        rsync -avzuP globus.opensciencedatacloud.org::public/illumina/igenomes/Homo_sapiens/UCSC/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa {DATA_DIR}
        mv {DATA_DIR}genome.fa {output}
        '''

rule index_genome:
	input: bwa = SRC_DIR + 'bwa/bwa',
               genome = DATA_DIR + 'hg19.fa'
	output: DATA_DIR + 'hg19.fa.bwt'
	params: h_vmem = '8g', bigio = '0',
            name = 'index_genome'
	log: LOG_DIR
	shell: '{input.bwa} index {input.genome}'

# Map with bwa (use backtrack algorithm because reads are 50 bp)
rule map_probes:
	input: fastq = DATA_DIR + 'ht12_probes.fq',
               genome = DATA_DIR + 'hg19.fa',
               bwa = SRC_DIR + 'bwa/bwa',
               index = DATA_DIR + 'hg19.fa.bwt'
	output: sai = DATA_DIR + 'ht12_probes.sai',
                sam = DATA_DIR + 'ht12_probes.sam'
	params: h_vmem = '12g', bigio = '0',
                name = 'map_probes'
	log: LOG_DIR
	shell: '''
        {input.bwa} aln {input.genome} {input.fastq} > {output.sai}
        {input.bwa} samse {input.genome} {output.sai} {input.fastq} > {output.sam}
        '''

rule sam_to_bam:
	input: DATA_DIR + 'ht12_probes.sam'              
	output: DATA_DIR + 'ht12_probes.bam'
	params: h_vmem = '8g', bigio = '0',
                name = 'sam_to_bam'
	log: LOG_DIR
	shell: 'samtools view -Sb {input} > {output}'

# This step removes unmapped probes, which cannot be disabled
rule bam_to_bed:
	input: DATA_DIR + 'ht12_probes.bam'
	output: DATA_DIR + 'ht12_probes.bed'
	params: h_vmem = '8g', bigio = '0',
                name = 'bam_to_bed'
	log: LOG_DIR
	shell: 'bedtools bamtobed -cigar -i {input} > {output}'

################################################################################
# Intersect probes and SNPs
################################################################################

rule intersect_bed:
	input: probes = DATA_DIR + 'ht12_probes.bed',
               snps  = DATA_DIR + 'snps_EUR_1KG.bed'
	output: DATA_DIR + 'ht12_probes_snps_ceu_hg19_af.bed'
	params: h_vmem = '16g', bigio = '0',
                name = 'intersect_bed'
	log: LOG_DIR
	shell: 'bedtools intersect -loj -a {input.probes} -b {input.snps} > {output}'

# After using bedtools intersect with the -loj flag, each probe is
# included in the new file at least once. If there is no SNP in the
# probe, the columns are filled with periods for character columns and
# -1 for numeric columns. However, if a probe has more than one SNP,
# it will be listed on multiple lines. The ultimate goal is to remove
# probes which contain SNPs at a certain minor allele frequency
# (MAF). Thus for proper filtering, we need to only keep the SNP with
# the highest MAF for each probe.
rule reduce_probes:
	input: DATA_DIR + 'ht12_probes_snps_ceu_hg19_af.bed'
	output: DATA_DIR + 'ht12_probes_snps_ceu_hg19_af_reduced.bed'
	params: h_vmem = '8g', bigio = '0',
                name = 'reduce_probes'
	log: LOG_DIR
	run:
          probes = open(input[0], 'r')
          good = open(output[0], 'w')

          d = {}
          for line in probes:
              cols = line.strip().split('\t')
              probe_id = cols[3]
              if probe_id not in d.keys():
                  d[cols[3]] = [line]
              else:
                  d[cols[3]].append(line)

          for key in d:
              if len(d[key]) == 1:
                  good.write(d[key][0])
              else: # Choose SNP with highest MAF for a given probe
                  top_maf = 0
                  result = ''
                  for snp in d[key]:
                      snp_af = snp.split('\t')[11]
                      if snp_af == 'NA':
                          snp_af = 0
                      snp_af = float(snp_af)
                      if snp_af > 0.5:
                          snp_af = 1 - snp_af
                      if snp_af >= top_maf:
                          top_maf = snp_af
                          result = snp
                  good.write(result)

          probes.close()
          good.close()

rule filter_probes:
	input: probes = DATA_DIR + 'ht12_probes_snps_ceu_hg19_af_reduced.bed',
               manifest = DATA_DIR + 'HumanHT-12_V4_0_R2_15002873_B.txt'
	output: probes = DATA_DIR + 'ht12_probes_snps_ceu_hg19_af_' + str(MAF) + '_map_' + str(MAP_SCORE) + '.txt',
                problem = DATA_DIR + 'problem_probes.txt'
	params: h_vmem = '12g', bigio = '0',
                name = 'filter_probes'
	log: LOG_DIR
	shell: 'Rscript analyze_probes.R {input.probes} {input.manifest} {MAF} {MAP_SCORE} {output.probes} -p {output.problem}'



'''
To submit:
nohup snakemake -kps make.py -j 96 --ri -c "qsub -l h_vmem={params.h_vmem} -l bigio={params.bigio} -N {params.name} -V -j y -cwd -o {log}" &
'''


# Paths
DATA_DIR = 'data/'
SRC_DIR = 'software/'
LOG_DIR = 'log/'

# Settings
CHROM = [str(x) for x in range(1, 23)] + ['X', 'Y', 'M']

# Target rules
localrules: all

rule all:
	input: DATA_DIR + 'snps_ceu_hg19_af.bed',
               DATA_DIR + 'ht12_probes.bed'

# Workflow
rule setup:
	output: DATA_DIR, SRC_DIR, LOG_DIR
	shell: 'mkdir -p {output}'

rule download_genos:
	output: DATA_DIR + 'genotypes_chr{CHR}_CEU_r28_nr.b36_fwd.txt.gz'
	shell: '''
    wget http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/2010-08_phaseII+III/forward/{output}
    mv {output} {DATA_DIR}
    '''

rule geno_to_bed:
    input: DATA_DIR + 'genotypes_chr{CHR}_CEU_r28_nr.b36_fwd.txt.gz'
	output: DATA_DIR + 'genotypes_chr{CHR}_CEU_r28_nr.b36_fwd.bed'
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

    # Read in genotypes and write out in bed format
    for g in geno:
        g = str(g, encoding='utf8')
        g_cols = g.strip().split(' ')
        chrom = g_cols[2]
        start = str(int(g_cols[3]) - 1)
        end = g_cols[3]
        name = g_cols[0]
        score = '0'
        strand = '+'
        bed.write(chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + score + '\t' + strand + '\n')

    geno.close()
    bed.close()
    
rule combine_genos:
	input: expand(DATA_DIR + 'genotypes_chr{CHR}_CEU_r28_nr.b36_fwd.bed', CHR = CHROM)
	output: DATA_DIR + 'snps_ceu_hg18.bed'
	params: h_vmem = '8g', bigio = '0',
            name = 'combine_genos'
	log: LOG_DIR
	shell: 'cat {input} > {output}'

rule download_lift_over:
	output: exec = SRC_DIR + 'liftOver',
            chain = DATA_DIR + 'hg18ToHg19.over.chain.gz'
	params: h_vmem = '8g', bigio = '0',
            name = 'download_lift_over'
	log: LOG_DIR
	shell: '''
    wget -O {output.exec} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
    chmod +x {output.exec}
    wget -O {output.chain} http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
    '''

rule lift_over_genos:
	input: genos = DATA_DIR + 'snps_ceu_hg18.bed',
           chain = DATA_DIR + 'hg18ToHg19.over.chain.gz'
	output: genos = DATA_DIR + 'snps_ceu_hg19.bed',
            lost = DATA_DIR + 'snps_lost.txt'
	params: h_vmem = '8g', bigio = '0',
            name = 'lift_over_genos'
	log: LOG_DIR
	shell: '{SRC_DIR}liftOver {input.genos} {input.chain} {output.genos} {output.lost}'

rule download_freqs:
	output: DATA_DIR + 'allele_freqs_chr{CHR}_CEU_r28_nr.b36_fwd.txt.gz'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'download_freqs' + wildcards.CHR
	log: LOG_DIR
	shell: 'wget -O {output} http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III/allele_freqs_chr{wildcards.CHR}_CEU_r28_nr.b36_fwd.txt.gz'

rule extract_freqs:
    input: DATA_DIR + 'allele_freqs_chr{CHR}_CEU_r28_nr.b36_fwd.txt.gz'
	output: DATA_DIR + 'allele_freqs_ceu_{CHR}.txt'
	params: h_vmem = '8g', bigio = '0',
            name = lambda wildcards: 'extract_freqs' + wildcards.CHR
	log: LOG_DIR
	run:
    import gzip

    # Open connection and skip header line
    freq = gzip.open(input[0], 'r')
    freq.readline()

    # Open connection to output file
    out = open(output[0], 'w')
    
    # Extract allele frequencies
    for f in freq:
        f = str(f, encoding='utf8')
        f_cols = f.strip().split(' ')
        out.write(f_cols[0] + '\t' + f_cols[11] + '\n')

    freq.close()
    out.close()

rule combine_freqs:
    input: expand(DATA_DIR + 'allele_freqs_ceu_{CHR}.txt', CHR = CHROM)
	output: DATA_DIR + 'allele_freqs_ceu.txt'
	params: h_vmem = '8g', bigio = '0',
            name = 'combine_freqs'
	log: LOG_DIR
	shell: 'cat {input} > {output}'

rule add_freqs:
	input: genos = DATA_DIR + 'snps_ceu_hg19.bed',
           freqs = DATA_DIR + 'allele_freqs_ceu.txt'
	output: DATA_DIR + 'snps_ceu_hg19_af.bed'
	params: h_vmem = '8g', bigio = '0',
                name = 'add_freqs'
	log: LOG_DIR
	run:
          freq = open(input.freqs, 'r')
          geno = open(input.genos, 'r')
    
          # Open connection to output file
          out = open(output[0], 'w')
    
          # Read in allele frequencies
          d_freq = {}
          for f in freq:
              f_cols = f.strip().split('\t')
              d_freq[f_cols[0]] = f_cols[1]

          for g in geno:
              g_cols = g.strip().split('\t')
              af = d_freq.get(g_cols[3], 'NA')
              out.write('\t'.join(g_cols[:4]) + '\t' + af + '\t' + g_cols[4] + '\n')

          freq.close()
          geno.close()
          out.close()

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
              if cols[0] == 'Homo sapiens':
                  name = '@' + cols[13] + ':' + cols[2] + ':' + cols[11]
                  seq = cols[17]
                  qual = '~' * len(seq)
                  fastq.write(name + '\n' + seq + '\n' + '+\n' + qual + '\n')

          manifest.close()
          fastq.close()

# Map with bwa (use backtrack algorithm because reads are 50 bp)
rule map_probes:
	input: fastq = DATA_DIR + 'ht12_probes.fq',
               genome = DATA_DIR + 'hg19.fa'
	output: sai = DATA_DIR + 'ht12_probes.sai',
                sam = DATA_DIR + 'ht12_probes.sam'
	params: h_vmem = '12g', bigio = '0',
                name = 'map_probes'
	log: LOG_DIR
	shell: '''
        software/bwa/bwa aln {input.genome} {input.fastq} > {output.sai}
        software/bwa/bwa samse {input.genome} {output.sai} {input.fastq} > {output.sam}
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


#probes_file <- commandArgs(trailingOnly = TRUE)[1]
probes_file <- "data/ht12_probes_snps_ceu_hg19_af_reduced.bed"
all_probes <- read.table(probes_file, sep = "\t", stringsAsFactors = FALSE)
colnames(all_probes) <- c("chr", "start", "end", "probe", "map_score", "strand",
                          "cigar", "chr_snp", "start_snp", "end_snp", "snp",
                          "snp_af", "snp_score")
# Check that each probe is only present in the file once
stopifnot(nrow(all_probes) == length(unique(all_probes$probe)))

################################################################################
# Plot dependence on mapping quality
################################################################################

# plot(density(all_probes$map_score))
# plot(all_probes$map_score)
# summary(all_probes$map_score)
qual_score_cutoff <- seq(min(all_probes$map_score),
                         max(all_probes$map_score),
                         by = 1)
num_probes_w_qual <- numeric(length = length(qual_score_cutoff))
for (i in seq_along(qual_score_cutoff)) {
  num_probes_w_qual[i] <- sum(all_probes$map_score >= qual_score_cutoff[i])
}
plot(qual_score_cutoff, num_probes_w_qual,
     xlab = "Mapping quality score cutoff",
     ylab = "Number of probes",
     main = "Probes vs. mapping quality score",
     type = "b")

################################################################################
# Plot dependence on SNP minor allele frequency (MAF)
################################################################################

af <- all_probes$snp_af
af <- ifelse(af > 0.5, 1 - af, af)

# sum(is.na(af)) # 14 probes have SNPs with no allele frequency estimate.
                 # These are not counted.
# Probes with no SNP have value of -1, thus they will always be lower than
# the minor allele frequency cutoff.

maf_cutoff <- seq(0, 0.5, by = 0.01)
num_probes_wo_snp <- numeric(length = length(maf_cutoff))
for (i in seq_along(maf_cutoff)) {
  num_probes_wo_snp[i] <- sum(af <= maf_cutoff[i], na.rm = TRUE)
}

plot(maf_cutoff, num_probes_wo_snp,
     xlab = "MAF cutoff",
     ylab = "Number of probes",
     main = "Probes vs. minor allele frequency",
     type = "b")

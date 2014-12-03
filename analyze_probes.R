suppressPackageStartupMessages(library(argparse))

get_args <- function(a) {
  ## Run 
  ## ./analyze_probes.R -h 
  ## for explanation of command line options
  parser <- ArgumentParser(description  = "Analyze and filter probes.")
  parser$add_argument("probes", nargs = 1,
                      help = "Bed file with all mapped probes with their intersected SNPs.")
  parser$add_argument("bam", nargs = 1,
                      help = "Bam file with mapping results")
  parser$add_argument("maf", nargs = 1,
                      help = "SNPs below this minor allele frequency cutoff will be ignored when filtering probes")
  parser$add_argument("mapping", nargs = 1,
                      help = "Probes with mapping score lower than this cutoff will be removed.")
  return(parser$parse_args(a))
}

main <- function(args) {
  probes_file <- args$probes
  all_probes <- read.table(probes_file, sep = "\t", stringsAsFactors = FALSE)
  colnames(all_probes) <- c("chr", "start", "end", "probe", "map_score", "strand",
                            "cigar", "chr_snp", "start_snp", "end_snp", "snp",
                            "snp_af", "snp_score")
  # Check that each probe is only present in the file once
  stopifnot(nrow(all_probes) == length(unique(all_probes$probe)))
  
  plot_probes_v_mapping(all_probes)
  plot_probes_v_maf(all_probes)
  
  # plot numbers lost
  # filter
  # add ensembl genes names
  # output
}


################################################################################
# Plot dependence on mapping quality
################################################################################

plot_probes_v_mapping <- function(all_probes) {
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
  png("probes_v_mapping_quality.png")
  plot(qual_score_cutoff, num_probes_w_qual,
       xlab = "Mapping quality score cutoff",
       ylab = "Number of probes",
       main = "Probes vs. mapping quality score",
       type = "b")
  dev.off()
}

################################################################################
# Plot dependence on SNP minor allele frequency (MAF)
################################################################################

plot_probes_v_maf <- function(all_probes) {
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
  png("probes_v_maf.png")
  plot(maf_cutoff, num_probes_wo_snp,
       xlab = "MAF cutoff",
       ylab = "Number of probes",
       main = "Probes vs. minor allele frequency",
       type = "b")
  dev.off()
}


if (!interactive()) {
  cline <- commandArgs(trailingOnly = TRUE)
} else {
  cline <- c("data/ht12_probes_snps_ceu_hg19_af_reduced.bed",
             "data/ht12_probes.bam",
             "0.05", "37")
}

args <- get_args(cline)
main(args)

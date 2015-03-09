suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))

get_args <- function(a) {
  ## Run 
  ## ./analyze_probes.R -h 
  ## for explanation of command line options
  parser <- ArgumentParser(description  = "Analyze and filter probes.")
  parser$add_argument("probes", nargs = 1,
                      help = "Bed file with all mapped probes with their intersected SNPs")
  parser$add_argument("manifest", nargs = 1,
                      help = "Illumina array manifest")
  parser$add_argument("maf", nargs = 1, type = "double",
                      help = "SNPs below this minor allele frequency cutoff will be ignored when filtering probes")
  parser$add_argument("mapping", nargs = 1, type = "integer",
                      help = "Probes with mapping score lower than this cutoff will be removed")
  parser$add_argument("output", nargs = 1, type = "character",
                      help = "File to write probe list")
  parser$add_argument("-p", "--problem", nargs = 1, type = "character", metavar = "filename",
                      help = "Write probes with ambiguous association with Ensembl ID to this file")
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
  
  # Calculate minor allele frequency
  all_probes$snp_maf <- ifelse(all_probes$snp_af > 0.5,
                               1 - all_probes$snp_af, all_probes$snp_af)
  
  test_snp_probe(all_probes)
  plot_probes_v_mapping(all_probes)
  plot_probes_v_maf(all_probes)
  filtered_probes <- filter_probes(all_probes, args$maf, args$mapping)
  formatted_probes <- format_probes(filtered_probes, args$problem)
  count_probes(all_probes, args$manifest, args$maf, args$mapping,
               formatted_probes)
  sorted_probes <- formatted_probes[order(formatted_probes$chr,
                                          formatted_probes$start,
                                          formatted_probes$end), ]
  write.table(sorted_probes, args$output,
              sep = "\t", quote = FALSE, row.names = FALSE)
}

################################################################################
# Test SNP-probe combinations
################################################################################

test_snp_probe <- function(all_probes) {
  ids  <- strsplit(all_probes$probe, ":")
  probeID <- sapply(ids, function(x) x[1])
  
  test_probe <- c("ILMN_2047240", "ILMN_2310685", "ILMN_2132898", "ILMN_2404154")
  test_snp <- c("rs6151429", "rs3794713", "rs1329151", "rs1303")
  test_chr <- c("chr22", "chr17", "chr10", "chr14")
  test_pos <- c(51063477, 80561616, 135234393, 94844843)
  for (i in seq_along(test_probe)) {
    if(all_probes[probeID == test_probe[i], "chr_snp"] != test_chr[i] &
       all_probes[probeID == test_probe[i], "end_snp"] != test_pos[i]) {
      warning(sprintf("Check versions of databases used.
                      Expect SNP %s (%s: %d) to fall in probe %s",
                      test_snp[i], test_chr[i], test_pos[i], test_probe[i]))
    }
  }
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
  
  # sum(is.na(af)) # 14 probes have SNPs with no allele frequency estimate.
  # These are not counted.
  # Probes with no SNP have value of -1, thus they will always be lower than
  # the minor allele frequency cutoff.
  
  maf_cutoff <- seq(0, 0.5, by = 0.01)
  num_probes_wo_snp <- numeric(length = length(maf_cutoff))
  for (i in seq_along(maf_cutoff)) {
    num_probes_wo_snp[i] <- sum(all_probes$snp_maf <= maf_cutoff[i],
                                na.rm = TRUE)
  }
  png("probes_v_maf.png")
  plot(maf_cutoff, num_probes_wo_snp,
       xlab = "MAF cutoff",
       ylab = "Number of probes",
       main = "Probes vs. minor allele frequency",
       type = "b")
  dev.off()
}

################################################################################
# Count number of probes after given filters
################################################################################

count_probes <- function(all_probes, manifest, maf, mapping, ens_probes) {

  command <- paste("grep 'Homo sapiens'", manifest, "| wc -l")
  total_num <- system(command, intern = TRUE)
  cat("Number of probes designed for human genome:", total_num, "\n")
  cat("Number of probes mapped to hg19:", nrow(all_probes), "\n")
  cat("Number of probes mapped uniquely to hg19:",
      sum(all_probes$map_score > 0), "\n")
  cat("Number of probes mapped with quality score >=", mapping, "to hg19:",
      sum(all_probes$map_score >= mapping), "\n")
  cat("Number of probes mapped without a SNP ( MAF >=", maf, "):",
      sum(all_probes$snp_maf <= maf, na.rm = TRUE), "\n")
  cat("Number of probes mapped with quality score >=", mapping,
      "to hg19 and without a SNP ( MAF >=", maf, "):",
      sum(all_probes$map_score >= mapping & all_probes$snp_maf <= maf,
          na.rm = TRUE), "\n")
  cat("Number of probes associated with a unique Ensembl gene ID",
      nrow(ens_probes), "\n")
}

################################################################################
# Filter probes
################################################################################

filter_probes <- function(all_probes, maf, mapping) {
  # Need to use function which when filtering to properly remove probes
  # with SNPs with allele frequency of NA
  filtered_probes <- all_probes[which(all_probes$map_score >= mapping &
                                      all_probes$snp_maf <= maf), ]
  return(filtered_probes)
}

################################################################################
# Format probes
################################################################################

format_probes <- function(probes, problem_file = NULL) {
  # problem_file - filename to save problem probes that are associated with more than
  # one Ensembl gene ID
  ensembl <- useMart(host = "feb2014.archive.ensembl.org",
                     biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl")
  ids  <- strsplit(probes$probe, ":")
  probeID <- sapply(ids, function(x) x[1])
  refseq <- sapply(ids, function(x) x[2])
  entrez <- sapply(ids, function(x) x[3])
  hugo <- sapply(ids, function(x) x[4])
  stopifnot(length(probeID) + length(refseq) + length(entrez) ==
            length(hugo) * 3)
  
  ens_gene_names <- getBM(attributes = c("ensembl_gene_id",
                                         "illumina_humanht_12_v4",
                                         "source", "chromosome_name",
                                         "hgnc_symbol",
                                         "gene_biotype",
                                         "entrezgene",
                                         "external_gene_id"),
                          filters = c("illumina_humanht_12_v4",
                                      "chromosome_name", "source"),
                          values = list(probeID, c(1:22, "X", "Y", "M"),
                                        c("ensembl", "ensembl_havana", "havana")),
                          mart = ensembl)
  
  ensg <- vector(length = nrow(probes))
  biotype <- vector(length = nrow(probes))
  ens_name <- vector(length = nrow(probes))

  problem <- vector()
  for (i in seq_along(probeID)) {
    probe_match <- which(ens_gene_names$illumina_humanht_12_v4 == probeID[i])
    if (entrez[i] != "") {
      entrez_match <- which(ens_gene_names$entrezgene == entrez[i])
    } else {
      entrez_match <- 1:nrow(probes)
    }
    if (!is.na(hugo[i])) {
      hugo_match <- which(ens_gene_names$hgnc_symbol == hugo[i])
    } else {
      hugo_match <- 1:nrow(probes)
    }
    row_match <- intersect(probe_match, intersect(entrez_match, hugo_match))
    result <- ens_gene_names[row_match, ]
    
    if (nrow(result) == 0) {
      ensg[i] <- NA
      biotype[i] <- NA
      ens_name[i] <- NA
    } else if (nrow(result) == 1) {
      ensg[i] <- result$ensembl_gene_id
      biotype[i] <- result$gene_biotype
      ens_name[i] <- result$external_gene_id
    } else if (length(unique(result$ensembl_gene_id)) == 1) {
      ensg[i] <- result$ensembl_gene_id[1]
      biotype[i] <- result$gene_biotype[1]
      ens_name[i] <- result$external_gene_id[1]
    } else {
      problem <- c(problem, i)
    }
  }
  
  if (!is.null(problem_file)) {
    df_problem <- vector()
    for (i in problem) {
      probe_match <- which(ens_gene_names$illumina_humanht_12_v4 == probeID[i])
      if (entrez[i] != "") {
        entrez_match <- which(ens_gene_names$entrezgene == entrez[i])
      } else {
        entrez_match <- 1:nrow(probes)
      }
      if (!is.na(hugo[i])) {
        hugo_match <- which(ens_gene_names$hgnc_symbol == hugo[i])
      } else {
        hugo_match <- 1:nrow(probes)
      }
      row_match <- intersect(probe_match, intersect(entrez_match, hugo_match))
      result <- ens_gene_names[row_match, ]
      df_problem <- rbind(df_problem, result)
    }
    write.table(df_problem, problem_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
  }
  final <- data.frame(chr = probes$chr,
                      start = probes$start,
                      end = probes$end,
                      probeID,
                      score = probes$map_score,
                      strand = probes$strand,
                      targetENSG = ensg,
                      targetHGNC = hugo,
                      targetEntrez = entrez,
                      targetBiotype = biotype,
                      targetENSName = ens_name)
  # Remove probes which are not associated with an Ensembl gene ID
  final <- final[!is.na(final$targetENSG), ]
  # Remove probes with ambiguous association with more than one Ensembl gene ID
  final <- final[final$targetENSG != FALSE, ]
  return(final)
}

if (!interactive()) {
  cline <- commandArgs(trailingOnly = TRUE)
} else {
  cline <- c("data/ht12_probes_snps_EUR_1KG_reduced_sorted.bed",
             "data/HumanHT-12_V4_0_R2_15002873_B.txt",
             "0.05", "37", "data/answer.txt")
}

args <- get_args(cline)
main(args)

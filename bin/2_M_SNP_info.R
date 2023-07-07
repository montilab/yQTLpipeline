#!/usr/bin/env Rscript

## calculate freq from gds file ##
## input description:
##  1 2. gds_name gds_file  name and file path for gds
##  3. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R
##  4. snp_assoc_txtfile  file contains SNP ids to test association; use "NA" to use all SNPs

## output:
## SNPinfo_{gds_name}.csv

args <- commandArgs(trailingOnly = TRUE)

gds_name <- args[1]
gds_file <- args[2]
shared_sampleid <- as.character(readLines(args[3]))
snp_assoc_txtfile <- args[4]

sink(paste0("2_", gds_name, "_freq.log"), append = FALSE, split = TRUE)
date()

suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(dplyr))

gds <- seqOpen(gds_file)

if (snp_assoc_txtfile != "NA") {
  snpid_df <- data.frame(
    "variant.id" = seqGetData(gds, "variant.id"),
    "annot.id" = seqGetData(gds, "annotation/id")
  )
  snpid_subset_annotid <- as.character(readLines(snp_assoc_txtfile))
  snpid_subset_varid <- snpid_df %>%
    filter(annot.id %in% snpid_subset_annotid) %>%
    pull(variant.id)
} else {
  snpid_subset_varid <- seqGetData(gds, "variant.id")
}

if (length(snpid_subset_varid) > 0) {

  seqSetFilter(gds, sample.id = shared_sampleid, variant.id = snpid_subset_varid)

  #### SNP allele info =========================================================
  SNP_info <- data.frame(
    "variant.id" = seqGetData(gds, "variant.id"),
    "REF" = seqGetData(gds, "$ref"),
    "ALT.ca" = seqGetData(gds, "$alt"),
    "chr" = seqGetData(gds, "chromosome"),
    "pos" = seqGetData(gds, "position"),
    "snpID" = seqGetData(gds, "annotation/id")
  )

  #### SNP miss + freq + MAC info ==============================================

  # order snpid_subset_varid so the allele frequency below matches the variant id order
  snpid_subset_varid <- snpid_subset_varid[order(snpid_subset_varid)]

  freq_df <- data.frame(
    "variant.id" = snpid_subset_varid,
    "freq.ca" = seqAlleleFreq(gds, ref.allele = 1L, minor = FALSE, .progress = FALSE, parallel = seqGetParallel(), verbose = FALSE),
    "freq.minor" = seqAlleleFreq(gds, ref.allele = 1L, minor = TRUE, .progress = FALSE, parallel = seqGetParallel(), verbose = FALSE),
    "miss" = seqMissing(gds, .progress = FALSE, parallel = seqGetParallel(), verbose = FALSE)
  ) %>%
    mutate(
      n.obs = round(length(shared_sampleid) * (1 - miss)),
      MAC = round(n.obs * 2 * freq.minor)
    )

  SNP_info <- merge(
    x = SNP_info,
    y = freq_df, by = "variant.id", all = F
  )

  #### Write ===================================================================
  seqResetFilter(gds)

  write.csv(SNP_info, paste0("SNPinfo_", gds_name, ".csv"), row.names = F)
  
} else {

  SNPinfo_expected_colnames <- c("variant.id", "REF", "ALT.ca", "chr", "pos", "snpID", "freq.ca", "freq.minor", "miss", "n.obs", "MAC")
  SNP_info <- data.frame(matrix(ncol = length(SNPinfo_expected_colnames), nrow = 0))
  colnames(SNP_info) <- SNPinfo_expected_colnames
  write.csv(SNP_info, paste0("SNPinfo_noSNP_", gds_name, ".csv"), row.names = F)
  
}

closefn.gds(gds)
date()
sink()

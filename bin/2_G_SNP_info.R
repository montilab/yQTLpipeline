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
suppressPackageStartupMessages(library(SNPRelate))
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

  #### SNP allele info ===========================================================
  SNP_info <- data.frame(
    "variant.id" = seqGetData(gds, "variant.id"),
    "REF" = seqGetData(gds, "$ref"),
    "ALT.ca" = seqGetData(gds, "$alt"),
    "snpID" = seqGetData(gds, "annotation/id")
  )

  #### Write =====================================================================
  seqResetFilter(gds)

  write.csv(SNP_info, paste0("SNPinfo_", gds_name, ".csv"), row.names = F)
} else {
  SNPinfo_expected_colnames <- c("variant.id", "REF", "ALT.ca", "snpID")
  SNP_info <- data.frame(matrix(ncol = length(SNPinfo_expected_colnames), nrow = 0))
  colnames(SNP_info) <- SNPinfo_expected_colnames
  write.csv(SNP_info, paste0("SNPinfo_noSNP_", gds_name, ".csv"), row.names = F)
}

closefn.gds(gds)
date()
sink()

#!/usr/bin/env Rscript

## input:
## 1. vcf_name
## 2. vcf_filepath
## ourput: gds file

args <- commandArgs(trailingOnly = TRUE)

vcf_name <- args[1]
vcf_filepath <- args[2]

sink("0_vcf_to_gds.log", append = FALSE, split = TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))

gds_file_name <- paste0(vcf_name, ".gds")

seqVCF2GDS(vcf_filepath, gds_file_name, verbose = TRUE)

date()
sink()

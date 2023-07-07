#!/usr/bin/env Rscript

## input:
## 1. {params.outdir}
## 2. vcf_list
## output: gds_list.csv

args <- commandArgs(trailingOnly = TRUE)

sink("0_write_gds_list.log", append = FALSE, split = TRUE)
date()

out_dir <- paste0(args[1], "/1_data/")
vcf_names <- read.csv(args[2], header = FALSE)[, 1]

cat("## Obtained", length(vcf_names), "VCF files. \n")
cat("## Output directory:", out_dir, "  \n")

gds_list <- cbind(
  vcf_names,
  paste0(args[1], "/1_data/", vcf_names, ".gds")
)

write.table(gds_list, "gds_list.csv",
  sep = ",", append = TRUE,
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

date()
sink()

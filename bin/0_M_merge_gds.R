#!/usr/bin/env Rscript

## input: gds files to merge; or, csv file of gds list
## output: merged.gds

args <- commandArgs(trailingOnly = TRUE)

sink("0_merge_gds.log", append = FALSE, split = TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))

if (endsWith(tolower(args[1]), ".csv")) {
  cat("Input CSV file contains GDS file paths.\n")
  gds_list <- read.csv(args[1], header = FALSE)[, 2]
} else if (endsWith(tolower(args[1]), ".gds")) {
  cat("Input GDS file list found.\n")
  gds_list <- args
}

cat(
  length(gds_list), "GDS found. \n",
  "Start merging...\n\n"
)

seqMerge(gds.fn = gds_list, out.fn = "merged.gds")

cat("Finished. \n")

date()
sink()

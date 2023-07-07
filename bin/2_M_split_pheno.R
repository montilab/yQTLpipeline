#!/usr/bin/env Rscript

## split phenotype ##
## input description:
##  1. pheno_file  file path for phenotype data, either .csv, .rds, .txt; need to have a column "sample.id" to match with .gds
##  2. pheno_name_all  use "all" for all columns (except "sample.id") in pheno_file; or txt file path which contains all the phenotype names to run, one column; those phenotype names should match with the column names in pheno_file
##  3. max_phenotype  numeric, maximum number of phenotype to run in parallel; technically there's no limit, but running too much phenotypes in parallel takes a lot of memory and make results impossible to save
##  4. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R

## outputs:
## "phenodat_chunk{i}.rds" "1_split_phenotype.log"

sink("2_split_phenotype.log", append = FALSE, split = TRUE)
date()

args <- commandArgs(trailingOnly = TRUE)
pheno_file <- args[1]
pheno_name_all <- as.character(readLines(args[2]))
max_phenotype <- as.numeric(args[3])
shared_sampleid <- as.character(readLines(args[4]))

#### read phenodat ####
phenodat <- NULL
if (endsWith(pheno_file, ".txt")) {
  phenodat <- read.table(pheno_file, sep = "\t", header = TRUE, check.names = FALSE)
} else if (endsWith(pheno_file, ".csv")) {
  phenodat <- read.csv(pheno_file, check.names = FALSE)
} else if (endsWith(pheno_file, ".rds")) {
  phenodat <- readRDS(pheno_file)
} else {
  stop("phenotype data must be .txt tab delim, .csv, or .rds")
}
if (!"sample.id" %in% colnames(phenodat)) {
  stop("phenotype data must contain a column 'sample.id' to match with sample id in genotype data. ")
}

## subset phenodat by shared_sampleid and pheno_name_all
phenodat$sample.id <- as.character(phenodat$sample.id)
phenodat <- phenodat[which(phenodat$sample.id %in% shared_sampleid), ]

if (pheno_name_all == "all") {
  pheno_name_all <- colnames(phenodat)
  pheno_name_all <- pheno_name_all[-which(pheno_name_all == "sample.id")]
}

if (!all(pheno_name_all %in% colnames(phenodat))) {
  pheno_name_all <- intersect(pheno_name_all, colnames(phenodat))
  warning("Some phenotype names specified are not found in phenotype data file. \n")
}

#### split ####
if (length(pheno_name_all) <= max_phenotype) {
  cat(paste(
    length(pheno_name_all),
    "phenotypes found.\nMaximum number of phenotype run in parallel:",
    max_phenotype, ".\n",
    "No need to split phenotype. \n"
  ))
  phenodat <- phenodat[, c("sample.id", pheno_name_all)]
  saveRDS(phenodat, "phenodat_chunk1.rds")
} else {
  pheno_name_chunklist <- split(pheno_name_all, ceiling(seq_along(pheno_name_all) / max_phenotype))
  cat(paste(
    length(pheno_name_all),
    "phenotypes found.\nMaximum number of phenotype run in parallel:",
    max_phenotype, "\nNeed to split into",
    length(pheno_name_chunklist), "chunks. \n\n"
  ))
  for (i in c(1:length(pheno_name_chunklist))) {
    cat(paste("  start spliting pheno chunk", i, "...\n"))
    pheno_name0 <- pheno_name_chunklist[[i]]
    phenodat0 <- phenodat[, c("sample.id", pheno_name0)]
    saveRDS(phenodat0, paste0("phenodat_chunk", i, ".rds"))
    cat(paste("  pheno chunk", i, "saved.\n"))
  }
}

cat(paste("\nFinished.\n"))

date()
sink()

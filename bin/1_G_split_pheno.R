#!/usr/bin/env Rscript

## split phenotype ##
## input description:
##  1. pheno_file  file path for phenotype data, either .csv, .rds, .txt; need to have a column "sample.id" to match with .gds
##  2. pheno_name_all  use "all" for all columns (except "sample.id") in pheno_file; or txt file path which contains all the phenotype names to run, one column; those phenotype names should match with the column names in pheno_file
##  3. shared_sampleid_file  shared sample id file, txt file generated from 2_get_sharedsampleid.R

## outputs:
## "phenodat_{phenoname0}.rds" "1_split_phenotype.log"

sink("1_split_phenotype.log", append = FALSE, split = TRUE)
date()

args <- commandArgs(trailingOnly = TRUE)
pheno_file <- args[1]
if (args[2] == "all") {
  pheno_name_all <- colnames(phenodat)
  pheno_name_all <- pheno_name_all[-which(pheno_name_all == "sample.id")]
} else {
  pheno_name_all <- readLines(args[2])
}
shared_sampleid <- as.character(readLines(args[3]))

#### read phenodat ####
phenodat <- NULL
if (endsWith(pheno_file, ".txt")) {
  phenodat <- read.table(pheno_file, sep = "\t", header = TRUE)
} else if (endsWith(pheno_file, ".csv")) {
  phenodat <- read.csv(pheno_file)
} else if (endsWith(pheno_file, ".rds")) {
  phenodat <- readRDS(pheno_file)
} else {
  stop("phenotype data must be .txt tab delim, .csv, or .rds")
}
if (!"sample.id" %in% colnames(phenodat)) {
  stop("phenotype data must contain a column 'sample.id' to match with sample id in genotype data. ")
}

## subset phenodat by shared_sampleid
phenodat$sample.id <- as.character(phenodat$sample.id)
phenodat <- phenodat[which(phenodat$sample.id %in% shared_sampleid), ]

#### split ####
## prepare pheno_name_all
if (!all(pheno_name_all %in% colnames(phenodat))) {
  pheno_name_all <- intersect(pheno_name_all, colnames(phenodat))
  warning("Warning: Some phenotype names specified are not found in phenotype data file. \n
  Please check if input phenotype names match column names in phenotype data. \n
  Or, use 'ALL' to automatically extract all columns in phenotype data. \n")
}

## get each phenotype data
if (length(pheno_name_all) >= 1) {
  cat(paste("\n", length(pheno_name_all), "phenotypes found.\n\n
Start preparing phenotype data...\n\n"))

  for (i in c(1:length(pheno_name_all))) {
    phenoname0 <- pheno_name_all[i]
    phenodat0 <- phenodat[, c("sample.id", phenoname0)]
    phenodat0 <- unique(phenodat0)
    saveRDS(phenodat0, paste0("phenodat_", phenoname0, ".rds"))
    writeLines(phenoname0, paste0("phenoname_", phenoname0, ".txt"))
    cat(paste("  Phenotype", i, ":", phenoname0, "saved.\n"))
  }
  cat(paste("\nAll pohenotype data saved.\nFinish.\n"))
} else {
  cat("Error: No phenotype found. Please check input phenotype data. \n\n
    Finish.\n")
}

date()
sink()

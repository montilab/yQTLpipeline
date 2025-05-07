#!/usr/bin/env Rscript

## get shared sampleid ##
## note: not subseting geno pheno cvrt file here;
##       pheno and cvrt will be subset by shared_sampleid in the split_pheno and prepare_cvrt step;
##       genodat will be subset in QTL analysis step.

## input description:
##  1. gds_file  file path for .gds
##  2. grm_file  "NA" if no grm
##  3. pheno_file  file path for phenotype data, either .csv, .rds, .txt; need to have a column "sample.id" to match with .gds
##  4. cvrt_file  "NA" if no covariates; file path for covariates data, either .csv, .rds, .txt; need to have a column "sample.id" to match with .gds
##  5. sample_id  "NA" if no; user specified sample id file

## outputs:
## "shared_sampleid.txt"  "0_sampleid.log"

sink("0_sampleid.log")
date()

args <- commandArgs(trailingOnly = TRUE)
gds_file <- args[1]
grm_file <- args[2]
pheno_file <- args[3]
cvrt_file <- args[4]
userdef_sampleid_file <- args[5]

suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(SeqArray))

#### user define sample id, if any ####
if (userdef_sampleid_file != "NA") {
  userdef_sampleid <- readLines(userdef_sampleid_file)
} else {
  userdef_sampleid <- NA
}

#### read genodat ####
showfile.gds(closeall = TRUE)
# gds_file <- "/restricted/projectnb/necs/Vanessa_Analysis/SomaCent_2021/gds/ALLmerged.gds"
genodat <- seqOpen(gds_file)
genodat_sampleid <- seqGetData(genodat, "sample.id")
cat("-- Sample size in genotype data:", length(genodat_sampleid), "\n")

#### read grm ####
if (grm_file == "NA") {
  grm_sampleid <- genodat_sampleid
} else {
  grm_sampleid <- colnames(readRDS(grm_file))
  cat("-- Sample size in GRM:", length(grm_sampleid), "\n")
}

#### read phenodat ####
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

cat("-- Sample size in phenotype data:", length(phenodat$sample.id), "\n")

#### read cvrtdat ####
cvrtdat <- NULL
if (cvrt_file != "NA" && cvrt_file != pheno_file) {
  if (endsWith(cvrt_file, ".txt")) {
    cvrtdat <- read.table(cvrt_file, sep = "\t", header = TRUE)
  } else if (endsWith(cvrt_file, ".csv")) {
    cvrtdat <- read.csv(cvrt_file)
  } else if (endsWith(cvrt_file, ".rds")) {
    cvrtdat <- readRDS(cvrt_file)
  } else {
    stop("covariate data must be .txt tab delim, .csv, or .rds")
  }
  if (!"sample.id" %in% colnames(cvrtdat)) {
    stop("covariate data must contain a column 'sample.id' to match with sample id in genotype data. ")
  }
} else {
  cvrtdat <- phenodat
}

cat("-- Sample size in covariates data:", length(cvrtdat$sample.id), "\n")

#### get shared.sampleid ####
cat("\nStart intersecting samples...\n")
shared_sampleid <- intersect(genodat_sampleid, intersect(grm_sampleid, intersect(phenodat$sample.id, cvrtdat$sample.id)))
if (!is.na(userdef_sampleid)) {
  shared_sampleid <- intersect(shared_sampleid, userdef_sampleid)
}
if (length(shared_sampleid) == 0) {
  stop("Error: there are no shared sample id between input genotype, phenotype and covariate data. ")
} else {
  cat(paste(
    "Shared samples: \n",
    "There are", length(shared_sampleid), "shared samples. \n\n",
    "A few examples:", paste(shared_sampleid[1:5], collapse = ", "), "\n\n"
  ))
}

write(shared_sampleid, "shared_sampleid.txt", ncolumns = 1)

date()
sink()

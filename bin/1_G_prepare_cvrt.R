#!/usr/bin/env Rscript

## convert factor cvrt to numeric ##
## input description:
##  1. cvrt_file  "NA" if no covariates; file path for covariates data, either .csv, .rds, .txt; need to have a column "sample.id" to match with .gds
##  2. PCs rds "NA" if no PCs or PCs are in cvrt file
##  3. covariates   "NA" if no covariates; or string, multiple covariates separate by ","; remember to specify the PCs here besides other covariates you may have
##  4. covariates_factor   "NA" if none; or string, specify which covariates are factors; multiple factor covariates separate by ","
##  5. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R

## outputs:
## "1_cvrt.rds"

args <- commandArgs(trailingOnly = TRUE)

sink("1_covariates.log", append = FALSE, split = TRUE)
date()

if (args[1] != "NA" && args[3] != "NA") {
  cvrt_file <- args[1]
  input_covariates <- unlist(strsplit(as.character(args[3]), ",", fixed = TRUE))
  if (args[3] != "NA") {
    factor_covariates <- unlist(strsplit(as.character(args[4]), ",", fixed = TRUE))
  }
}

shared_sampleid <- as.character(readLines(args[5]))

## read cvrtdat
cvrtdat <- NULL
if (cvrt_file != "NA" && input_covariates != "NA") {
  if (endsWith(tolower(cvrt_file), ".txt")) {
    cvrtdat <- read.table(cvrt_file, sep = "\t", header = TRUE)
  } else if (endsWith(tolower(cvrt_file), ".csv")) {
    cvrtdat <- read.csv(cvrt_file)
  } else if (endsWith(tolower(cvrt_file), ".rds")) {
    cvrtdat <- readRDS(cvrt_file)
  } else {
    stop("Error: phenotype and covariates data must be in .txt/TXT tab delim, .csv/CSV, or .rds/RDS format.")
  }

  if (!"sample.id" %in% colnames(cvrtdat)) {
    stop("Error: phenotype and covariates data must contain a column 'sample.id' to match with the sample id in genotype data. ")
  }
}

if (args[2] != "NA") {
  PC_rds <- readRDS(args[2])
  duplicate_cols <- intersect(colnames(PC_rds), colnames(cvrtdat))
  duplicate_cols <- duplicate_cols[which(duplicate_cols != "sample.id")]
  if (length(duplicate_cols) > 0) {
    cat(
      "Warning: Duplicated columns found in input covariate file and PCs RDS file:\n",
      paste(duplicate_cols, collapse = ", "),
      "\nDiscard those columns in the covariates file. \n"
    )
    cvrtdat[, duplicate_cols] <- NULL
  }
  cvrtdat <- merge(x = cvrtdat, y = PC_rds, by = "sample.id", all = FALSE)
}

cvrtdat <- subset(cvrtdat, sample.id %in% shared_sampleid)
covariates <- intersect(colnames(cvrtdat), input_covariates)
factor_covariates <- intersect(covariates, factor_covariates)
numeric_covariates <- setdiff(covariates, factor_covariates)
cvrtdat <- cvrtdat[, c("sample.id", covariates)]

for (n in numeric_covariates) {
  cvrtdat[, n] <- as.numeric(cvrtdat[, n])
}
cat("\nSince the algorithm only takes numeric covariate, factor covariates with
    multiple levels will be converted into multiple columns.\n\n")

if (length(factor_covariates) >= 1) {
  cat(paste("Factor covariates:", paste(factor_covariates, collapse = ", "), ".\n"))
  for (factor_name0 in factor_covariates) {
    cat(paste("-- Start processing", factor_name0, "...\n"))
    f_data0 <- cvrtdat[, factor_name0]
    all_levels0 <- as.character(unique(f_data0))
    if (length(all_levels0) > 2) {
      cat(paste(
        "  ", factor_name0, "has", length(all_levels0), "levels. Including:\n ",
        paste(all_levels0, collapse = ", "),
        "\n    Set", all_levels0[1], "as the baseline.\n"
      ))
      f_data_mat0 <- matrix(
        nrow = nrow(cvrtdat),
        ncol = length(all_levels0) - 1,
        dimnames = list(cvrtdat$sample.id, paste(factor_name0, all_levels0[2:length(all_levels0)], sep = "_"))
      )
      for (l in all_levels0[2:length(all_levels0)]) {
        # set the first level as the "baseline"
        f_data_mat0[, paste(factor_name0, l, sep = "_")] <- ifelse(test = f_data0 == l, yes = 1, no = 0)
      }
      cvrtdat[, factor_name0] <- NULL
      cvrtdat <- cbind(cvrtdat, f_data_mat0)
      cat(paste("-- Finish processing", factor_name0, ".\n"))
    } else {
      cat(paste("  ", factor_name0, "only has 2 levels. \n"))
      cvrtdat[, factor_name0] <- as.numeric(as.factor(cvrtdat[, factor_name0])) - 1
      cat("   Converted into 0 and 1. \n")
    }
  }
  cat("\nAll factor covariates finished. \n")
} else {
  cat("\nNo factor covariate need processing. \n")
}

cat("Covariate used: ", paste(covariates, collapse = ", "))

saveRDS(cvrtdat, "1_cvrt.rds")
cat(paste("\nCovariates data saved.\n"))

date()
sink()

#!/usr/bin/env Rscript

## null model and association test ##
## input description:
##  1 2. gds_name gds_file  name and file path for gds
##  3. "genotype" for genotype data. "dosage" to use dosage (i.e. impute data)
##  4. cvrt_rds  rds file path for cvrt .rds; has column "sample.id" to match with .gds
##  5. grm_rds  rds file path for grm matrix
##  6. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R
##  7. snp_assoc_txtfile  file contains SNP ids to test association; use "NA" to use all SNPs
##  8. pheno_rds  rds file path for one single phenotype data; has column "sample.id" to match with .gds
##  9. pval_cutoff  p-value cutoff for QTL

## output:
## nullmod_{phenoname}.rds  assoc_{phenoname}.rds  count_{phenoname}.txt

args <- commandArgs(trailingOnly = TRUE)
phenoname <- gsub(".*phenodat_(.+).rds.*", "\\1", args[8])

gds_name <- args[1]
gds_file <- args[2]

sink(paste0("3_", phenoname, "_", gds_name, "_assoc.log"), append = FALSE, split = TRUE)
date()

stopifnot(
  "Error: Parameter 'genotype_or_dosage' must be one of: 'GT', 'genotype', 'DS', 'dosage'. " =
    args[3] %in% c("GT", "genotype", "DS", "dosage")
)
use_impute <- ifelse(args[3] %in% c("DS", "dosage"), TRUE, FALSE)

cvrt_rds <- readRDS(args[4])
if (args[5] == "NA") {
    grm <- NULL 
} else {
    grm <- readRDS(args[5])
}
shared_sampleid <- as.character(readLines(args[6]))
snp_assoc_txtfile <- args[7]
pheno_dat <- readRDS(args[8])
pval_cutoff <- as.numeric(args[9])

suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(GENESIS))

showfile.gds(closeall = TRUE, verbose = FALSE)
gds <- seqOpen(gds_file)
run_assoc <- FALSE
assoc_has_result <- FALSE

#### First check if required SNP are in GDS file ===============================
if (snp_assoc_txtfile == "NA") {
  snpid_subset_varid <- seqGetData(gds, "variant.id")
} else {
  snpid_subset_annotid <- as.character(readLines(snp_assoc_txtfile))
  snpid_df <- data.frame(
    "variant.id" = seqGetData(gds, "variant.id"),
    "annot.id" = seqGetData(gds, "annotation/id")
  )
  snpid_subset_varid <- snpid_df %>%
    filter(annot.id %in% snpid_subset_annotid) %>%
    pull(variant.id)
}
if (length(snpid_subset_varid) > 0) {
  run_assoc <- TRUE
}

#### Read data and run association =============================================
if (run_assoc) {
  assoc_has_result <- FALSE
  seqSetFilter(gds, variant.id = snpid_subset_varid)

  ##### Create annot: merge cvrt (which has PC) + phenotype --------------------
  covariates <- colnames(cvrt_rds)
  covariates <- covariates[!covariates %in% "sample.id"]
  annot <- merge(x = cvrt_rds, y = pheno_dat, by = "sample.id", all.x = FALSE, all.y = FALSE) %>%
    mutate(sample.id = as.character(sample.id))

  gds_sampleid <- data.frame("sample.id" = as.character(seqGetData(gds, "sample.id")))
  annot <- left_join(
    x = gds_sampleid, y = annot,
    by = "sample.id"
  )
  annot_df <- AnnotatedDataFrame(annot)

  ##### SeqVarData object ------------------------------------------------------
  seqData <- SeqVarData(gds, sampleData = annot_df)

  ##### Null model -------------------------------------------------------------
  cat(paste("\n-- Begin Nullmod of", gds_name, phenoname, ".\n"))
  nullmod <- try(fitNullModel(seqData,
    outcome = phenoname,
    covars = covariates,
    cov.mat = grm,
    family = "gaussian",
    verbose = TRUE,
    sample.id = shared_sampleid
  ))

  ##### GWAS -------------------------------------------------------------------
  if (class(nullmod) != "try-error") {
    saveRDS(nullmod, paste0("nullmod_", phenoname, "_", gds_name, ".rds"))
    cat("\n-- Nullmod finished. \n\n")

    cat(paste("\n-- Begin association of", gds_name, phenoname, ".\n"))
    iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)

    assoc <- assocTestSingle(iterator, imputed = use_impute, nullmod, test = "Score", verbose = TRUE)
    assoc <- assoc %>% filter(Score.pval <= pval_cutoff)
    cat("\n-- Association finished. \n\n")

    if (nrow(assoc) > 0) {
      assoc_has_result <- TRUE

      ## Add missingness, phenotype, format columns
      assoc <- assoc %>%
        mutate(.after = "n.obs", miss = 1 - (assoc$n.obs / length(shared_sampleid))) %>%
        mutate(.before = everything(), phenotype = phenoname) %>%
        mutate(
          chr = as.numeric(as.character(chr)),
          pos = as.numeric(as.character(pos))
        ) %>% 
        rename_with(~paste0(., ".ca"), .cols = c(allele.index, freq, Score, Score.Stat, Est))

      saveRDS(assoc, paste0("assoc_", phenoname, "_", gds_name, ".rds"))
      counttable <- c(phenoname, gds_name, nrow(assoc), "success")
      cat("-- Finish saving", paste0("assoc_", phenoname, "_", gds_name, ".rds"), ".\n\n")
    } else {
      cat(paste0("Message: ", phenoname, gds_name, " do not have QTLs passed threshold.\n"))
      counttable <- c(phenoname, gds_name, 0, "success_noQTL")
      assoc_has_result <- FALSE
    }
    seqResetFilter(iterator)
  } else {
    cat(paste0("Warning: ", phenoname, gds_name, " failed to generate nullmodel.\n"))
    counttable <- c(phenoname, gds_name, 0, "failed")
    assoc_has_result <- FALSE
  }
  seqResetFilter(gds)
} else {
  cat(paste("Message:", gds_name, "do not contain any input SNPs.\n"))
  counttable <- c(phenoname, gds_name, 0, "no_input_SNP")
  assoc_has_result <- FALSE
}

write(counttable, paste0("count_", phenoname, "_", gds_name, ".txt"), sep = "\t", ncolumns = 4)
# Still needs to produce expected "count_*.txt" even if no results.

#### Write psedo result if no results ==========================================
if (!assoc_has_result | !run_assoc) {
  assoc_expected_colnames <- c("phenotype", "variant.id", "chr", "pos", "allele.index.ca", "n.obs", "miss", "freq.ca", "MAC", "Score.ca", "Score.SE", "Score.Stat.ca", "Score.pval", "Est.ca", "Est.SE", "PVE")
  assoc_psedo_res <- data.frame(matrix(ncol = length(assoc_expected_colnames), nrow = 0))
  colnames(assoc_psedo_res) <- assoc_expected_colnames
  assoc_psedo_res$chr <- as.numeric(assoc_psedo_res$chr)
  assoc_psedo_res$pos <- as.numeric(assoc_psedo_res$pos)
  saveRDS(assoc_psedo_res, paste0("assoc_noQTL_", phenoname, "_", gds_name, ".rds"))
  # Still needs to produce expected rds "assoc_*.rds" even if no results.
}

closefn.gds(gds)
date()
sink()

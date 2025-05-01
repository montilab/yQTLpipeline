#!/usr/bin/env Rscript

## null model and association test ##
## input description:
##  1 2. gds_name gds_file  name and file path for gds
##  3. "genotype" for genotype data. "dosage" to use dosage (i.e. impute data)
##  4. cvrt_rds  rds file path for cvrt .rds; has column "sample.id" to match with .gds
##  5. grm_rds  rds file path for grm matrix
##  6. nullmod_family  "gaussian"/"linear" or "binomial"/"logistic"
##  7. test_method  "Score" or "Score.SPA"
##  8. pval_cutoff  p-value cutoff for QTL
##  9. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R
## 10. snp_assoc_txtfile  file contains SNP ids to test association; use "NA" to use all SNPs
## 11. pheno_rds  rds file path for one single phenotype data; has column "sample.id" to match with .gds
## 12. draw_genopheno_boxplot  "true" to draw genotype-phenotype boxplot for top SNP in each chromosome
## 13. boxplot_p_cutoff  the p-value cutoff to draw genotype-phenotype boxplot

## output:
## nullmod_{phenoname}.rds  assoc_{phenoname}.rds  count_{phenoname}.txt

args <- commandArgs(trailingOnly = TRUE)
phenoname <- gsub(".*phenodat_(.+).rds.*", "\\1", args[11])

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

args[6] <- tolower(args[6])
stopifnot(
  "Error: Parameter 'nullmod_family' must be one of: 'gaussian' (or 'linear') or 'binomial' (or 'logistic'). " =
    args[6] %in% c("gaussian", "linear", "binomial", "logistic")
)
if (args[6] %in% c("gaussian", "linear")) nullmod_family <- "gaussian"
if (args[6] %in% c("binomial", "logistic")) nullmod_family <- "binomial"

args[7] <- tolower(args[7])
stopifnot(
  "Error: Parameter 'test_method' must be one of: 'Score' or 'Score.SPA' (or 'SPA'). " =
    args[7] %in% c("score", "score.spa", "spa")
)
if (args[7] == "score") test_method <- "Score"
if (args[7] %in% c("score.spa", "spa")) test_method <- "Score.SPA"

pval_cutoff <- as.numeric(args[8])
shared_sampleid <- as.character(readLines(args[9]))
snp_assoc_txtfile <- args[10]
pheno_dat <- readRDS(args[11])

draw_genopheno_boxplot <- FALSE
boxplot_p_cutoff <- 1
if (args[12] == "true") {
  draw_genopheno_boxplot <- TRUE
  boxplot_p_cutoff <- as.numeric(args[13])
}

suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
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
    family = nullmod_family,
    verbose = TRUE,
    sample.id = shared_sampleid
  ))

  ##### GWAS -------------------------------------------------------------------
  if (class(nullmod) != "try-error") {
    saveRDS(nullmod, paste0("nullmod_", phenoname, "_", gds_name, ".rds"))
    cat("\n-- Nullmod finished. \n\n")

    cat(paste("\n-- Begin association of", gds_name, phenoname, ".\n"))
    iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)

    assoc <- assocTestSingle(iterator, imputed = use_impute, nullmod, test = test_method, verbose = TRUE)

    if (test_method == "Score") {
      assoc <- assoc %>% filter(Score.pval <= pval_cutoff)
    } else if (test_method == "Score.SPA") {
      assoc <- assoc %>% filter(SPA.pval <= pval_cutoff)
    }
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
        rename_with(~ paste0(., ".ca"), .cols = c(allele.index, freq, Score, Score.Stat, Est))

      saveRDS(assoc, paste0("assoc_", phenoname, "_", gds_name, ".rds"))
      counttable <- c(phenoname, gds_name, nrow(assoc), "success")
      cat("-- Finish saving", paste0("assoc_", phenoname, "_", gds_name, ".rds"), ".\n\n")

      ##### Plot genotype-phenotype plot ---------------------------------------
      ## assoc has been saved alreayd. rename the p value column to "pval" is only used in this script for plotting
      assoc <- assoc %>% rename(pval = any_of(c("Score.pval", "SPA.pval")))

      if (draw_genopheno_boxplot & min(assoc$pval) < boxplot_p_cutoff) {
        cat("-- Start drawing phenotype-genotype boxplot. \n")
        ## first, obtain all SNPs passed p threshold
        assoc_sub <- assoc %>%
          select(variant.id, chr, pos, pval) %>%
          filter(pval < boxplot_p_cutoff)

        cat("-- Start subset. \n")
        ## subset GDS to only selected variants
        seqSetFilter(gds, variant.id = assoc_sub$variant.id, sample.id = shared_sampleid)

        cat("-- Start SNP info. \n")
        ## obtain SNP information
        SNP_info_sub <- data.frame(
          "variant.id" = seqGetData(gds, "variant.id"),
          "REF" = seqGetData(gds, "$ref"),
          "ALT" = seqGetData(gds, "$alt")
        )

        ## GENESIS only process one phenotype at a time, no need to loop through phenotypes
        ## for each chromosome, only keep the top SNP
        ## sometimes multiple SNPs have the same p-value. keep the smallest pos
        merged_sub <- merge(x = SNP_info_sub, y = assoc_sub, by = "variant.id", all = T)
        merged_sub_filtered <- merged_sub %>%
          group_by(chr) %>%
          arrange(pos) %>%
          slice_min(order_by = pval) %>%
          distinct(chr, .keep_all = TRUE) %>%
          as.data.frame()
        cat("-- Number of top SNPs:", nrow(merged_sub_filtered), "\n")

        ## obtain genotype information for plotting
        geno_mat_sub <- snpgdsGetGeno(gds,
          snpfirstdim = TRUE, with.id = TRUE,
          snp.id = merged_sub_filtered$variant.id, sample.id = shared_sampleid
        )
        colnames(geno_mat_sub$genotype) <- geno_mat_sub$sample.id
        rownames(geno_mat_sub$genotype) <- geno_mat_sub$snp.id
        geno_mat_sub <- as.data.frame(t(geno_mat_sub$genotype)) %>%
          mutate(across(everything(), ~ factor(case_when(
            . == 0 ~ "Alt/Alt",
            . == 1 ~ "Ref/Alt",
            . == 2 ~ "Ref/Ref",
            TRUE ~ as.character(.)
          ), levels = c("Ref/Ref", "Ref/Alt", "Alt/Alt")))) %>%
          rownames_to_column(var = "sample.id")

        ## merge genotype and phenotype data
        plot_pheno_dat <- merge(
          x = geno_mat_sub,
          y = pheno_dat, # pheno_dat only has two columns: sample.id, phenotype (phenoname)
          by = "sample.id", all.x = F, all.y = F
        )

        ## plot for each SNP
        for (variant0 in merged_sub_filtered$variant.id) {
          snp_name0 <- paste0("chr", paste(merged_sub_filtered %>% filter(variant.id == variant0) %>% select(chr, pos, REF, ALT) %>% as.character(), collapse = "_"))
          snp_column0 <- as.character(variant0)
          png(paste0(phenoname, "_", snp_name0, ".png"), width = 600, height = 400, res = 120)
          print(
            ggplot(plot_pheno_dat, aes(x = !!as.name(snp_column0), y = !!as.name(phenoname), fill = !!as.name(snp_column0))) +
              geom_boxplot() +
              theme_minimal() +
              xlab(snp_name0) +
              theme(legend.position = "none") +
              scale_fill_manual(values = c("#E5A1C3", "#A1A1E5", "#A1E5C3")) +
              ggtitle(paste(phenoname, "x", snp_name0))
          )
          dev.off()
        }
      }
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
  if (test_method == "Score") {
    assoc_expected_colnames <- c("phenotype", "variant.id", "chr", "pos", "allele.index.ca", "n.obs", "miss", "freq.ca", "MAC", "Score.ca", "Score.SE", "Score.Stat.ca", "Score.pval", "Est.ca", "Est.SE", "PVE")
  } else if (test_method == "Score.SPA") {
    assoc_expected_colnames <- c("phenotype", "variant.id", "chr", "pos", "allele.index.ca", "n.obs", "miss", "freq.ca", "MAC", "Score.ca", "Score.SE", "Score.Stat.ca", "SPA.pval", "Est.ca", "Est.SE", "PVE", "SPA.converged")
  }
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

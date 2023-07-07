#!/usr/bin/env Rscript

## QTL analysis ## matrixeQTL ##

## input description:
##  1 2. gds_name gds_file  file path for .gds
##  3. cvrt_rds  rds file path for cvrt .rds; has column "sample.id" to match with .gds
##  4. pval_cutoff  p-value cutoff for QTL
##  5. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R
##  6. snp_assoc_txtfile  file contains SNP ids to test association; use "NA" to use all SNPs
##  7. pheno_rds  rds file path for phenotype data (each chunk); has column "sample.id" to match with .gds

## output:
## "2_QTL_[gds_name]_[chunk_name].log"
## "QTL_[gds_name]_[chunk_name].rds"
## QTL_count_[gds_name]_[chunk_name].txt"

args <- commandArgs(trailingOnly = TRUE)
gds_name <- args[1]
gds_file <- args[2]
chunk_name <- gsub(".*phenodat_(.+).rds.*", "\\1", args[7])

sink(paste0("3_QTL_", gds_name, "_", chunk_name, ".log"), append = FALSE, split = TRUE)

cvrt_dat <- readRDS(args[3])
pval_cutoff <- as.numeric(args[4])
shared_sampleid <- as.character(readLines(args[5]))
snp_assoc_txtfile <- args[6]
pheno_dat <- readRDS(args[7])
phenonames <- colnames(pheno_dat)[which(colnames(pheno_dat) != "sample.id")]

date()

suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(MatrixEQTL))

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

if (run_assoc) {

  #### convert GDS to matrix =====================================================
  cat("## Start converting genotype data...\n")
  genomat <- snpgdsGetGeno(gds, snpfirstdim = TRUE, with.id = TRUE, snp.id = snpid_subset_varid, sample.id = shared_sampleid)

  colnames(genomat$genotype) <- genomat$sample.id
  rownames(genomat$genotype) <- genomat$snp.id
  cat("## Finish convertion. \n\n")

  #### sort according to shared.sample.id ======================================
  geno_mat <- genomat$genotype[, match(shared_sampleid, genomat$sample.id)]
  saveRDS(head(geno_mat), "input_geno_mat_header.rds")

  cvrt_mat <- cvrt_dat[match(shared_sampleid, cvrt_dat$sample.id), ]
  cvrt_mat$sample.id <- NULL
  rownames(cvrt_mat) <- shared_sampleid
  saveRDS(cvrt_mat, "input_cvrt_pc.rds")
  cat("## Finish reading covariates. \n\n")

  #### obtain phenotype data ===================================================

  cat("## Start analyzing", chunk_name, "...\n")

  pheno_mat <- pheno_dat[match(shared_sampleid, pheno_dat$sample.id), ]
  rownames(pheno_mat) <- shared_sampleid
  pheno_mat$sample.id <- NULL
  saveRDS(pheno_mat, paste0("input_phenotype_", chunk_name, ".rds"))

  cat("## Start matrixeQTL engine...\n")

  #### matrixeQTL ####
  useModel <- modelLINEAR
  output_file_name <- tempfile()
  errorCovariance <- numeric()

  # load SNP info
  snps <- SlicedData$new()
  snps$CreateFromMatrix(geno_mat)

  # load covariates info
  cvrt <- SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(t(cvrt_mat)))

  # load phenotype info
  gene <- SlicedData$new()
  gene$CreateFromMatrix(as.matrix(t(pheno_mat)))

  me <- Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pval_cutoff,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )

  cat("## Finish matrixeQTL engine. Continue to format the result...\n")

  ##### QTL result -------------------------------------------------------------
  QTL_res_df <- me$all$eqtls

  if (nrow(QTL_res_df) > 0) {
    assoc_has_result <- TRUE
    cat("## Finish QTL analysis.\n\n")

    ## rename columns; flip the direction of the effect estimate so it will be the effect of ALT
    QTL_res_df <- QTL_res_df %>%
      rename(variant.id = snps, phenotype = gene) %>%
      mutate(statistic = (-1) * statistic, beta = (-1) * beta, FDR = NULL) %>% 
      rename(beta.ca = beta, statistic.ca = statistic)

    saveRDS(QTL_res_df, paste0("QTL_", gds_name, "_", chunk_name, ".rds"))

    ##### QTL count table ------------------------------------------------------
    QTLcount <- QTL_res_df %>%
      count(phenotype, name = "QTL_num") %>%
      mutate(genotype_dat = gds_name, .before = QTL_num) %>%
      mutate(description = "success", .after = everything())

    ## some phenotypes may not have QTLs passed threshold
    phenotypes_no_QTL <- setdiff(phenonames, unique(QTL_res_df$phenotype))
    if (length(phenotypes_no_QTL) > 0) {
      cat("Message: some phenotypes in", chunk_name, gds_name, "do not have QTLs passed threshold. See QTL_counts_summary.txt for details. \n")
      QTLcount2 <- data.frame(
        "phenotype" = phenotypes_no_QTL,
        "genotype_dat" = gds_name,
        "QTL_num" = 0,
        "description" = "noQTL"
      )
      QTLcount <- rbind(QTLcount, QTLcount2)
    }

    cat("## Finish saving results for", gds_name, "_", chunk_name, ".\n\n")
  } else {
    assoc_has_result <- FALSE
    cat(paste("Message: all phenotypes in", chunk_name, gds_name, "do not have QTLs passed threshold.\n"))
    QTLcount <- data.frame(
      "phenotype" = phenonames,
      "genotype_dat" = gds_name,
      "QTL_num" = 0,
      "description" = "noQTL"
    )
  }
} else {
  cat(paste("Message:", gds_name, "do not contain any input SNPs.\n"))
  QTLcount <- data.frame(
    "phenotype" = phenonames,
    "genotype_dat" = gds_name,
    "QTL_num" = 0,
    "description" = "no_input_SNP"
  )
}

if (!run_assoc | !assoc_has_result) {
  QTL_expected_colnames <- c("variant.id", "phenotype", "statistic.ca", "pvalue", "beta.ca")
  QTL_psedo_res <- data.frame(matrix(ncol = length(QTL_expected_colnames), nrow = 0))
  colnames(QTL_psedo_res) <- QTL_expected_colnames
  saveRDS(QTL_psedo_res, paste0("QTL_noQTL_", gds_name, "_", chunk_name, ".rds"))
}

write.table(QTLcount, paste0("QTL_count_", gds_name, "_", chunk_name, ".txt"),
  quote = F, sep = "\t", col.names = T, row.names = F
)

cat("## Finish. \n")

closefn.gds(gds)
date()
sink()

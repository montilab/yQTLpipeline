#!/usr/bin/env Rscript

## QTL analysis ## matrixeQTL ##

## input description:
##  1 2. gds_name gds_file  file name and path for .gds
##  3. cvrt_rds  rds file path for cvrt .rds; has column "sample.id" to match with .gds
##  4. model_type  character input, one of: "linear", "category", "interaction"
##  5. interaction_variable  the cvrt to test interaction when model_type is "interaction". Use NA as place holder otherwise.
##  6. pval_cutoff  p-value cutoff for QTL
##  7. shared_sampleid_file  shared sample id file, txt file generated from 0_get_sampleid.R
##  8. snp_assoc_txtfile  file contains SNP ids to test association; use "NA" to use all SNPs
##  9. pheno_rds  rds file path for phenotype data (each chunk); has column "sample.id" to match with .gds
##  10. draw_genopheno_boxplot  "true" to draw genotype-phenotype boxplot for top SNP in each chromosome
##  11. boxplot_p_cutoff  the p-value cutoff to draw genotype-phenotype boxplot

## output:
## "2_QTL_[gds_name]_[chunk_name].log"
## "QTL_[gds_name]_[chunk_name].rds"
## "QTL_count_[gds_name]_[chunk_name].txt"
## "[phenoname]_[snp_name].png" if draw genotyoe-phenotype boxplot

args <- commandArgs(trailingOnly = TRUE)
gds_name <- args[1]
gds_file <- args[2]
chunk_name <- gsub(".*phenodat_(.+).rds.*", "\\1", args[9])

sink(paste0("3_QTL_", gds_name, "_", chunk_name, ".log"), append = FALSE, split = TRUE)
date()

suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MatrixEQTL))

cvrt_dat <- readRDS(args[3])
model_type <- args[4]
if (model_type == "linear") {
  useModel <- modelLINEAR
} else if (model_type == "category") {
  useModel <- modelANOVA
} else if (model_type == "interaction") {
  useModel <- modelLINEAR_CROSS
  interaction_cvrt <- args[5]
  if (is.na(interaction_cvrt) | interaction_cvrt == "NA") {
    stop("Error: model_type is specified as \"interaction\" but no interaction covariate specified. \n")
  }
  if (!interaction_cvrt %in% colnames(cvrt_dat)) {
    stop("Error: interaction covariate does not exist in the input covariate data. \n")
  }
  ## relocate the interaction covariate to the last column as required by using modelLINEAR_CROSS
  cvrt_dat <- cvrt_dat %>% relocate(!!as.name(interaction_cvrt), .after = everything())
} else {
  stop("Error: Parameter \"model_type\" must be one of the following: \"linear\", \"category\", \"interaction\". \n")
}
pval_cutoff <- as.numeric(args[6])
shared_sampleid <- as.character(readLines(args[7]))
snp_assoc_txtfile <- args[8]
pheno_dat <- readRDS(args[9])
phenonames <- colnames(pheno_dat)[which(colnames(pheno_dat) != "sample.id")]

if (args[10] == "true") {
  draw_genopheno_boxplot <- TRUE
  boxplot_p_cutoff <- as.numeric(args[11])
}

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
  cat("-- Start converting genotype data...\n")
  genomat <- snpgdsGetGeno(gds,
    snpfirstdim = TRUE, with.id = TRUE,
    snp.id = snpid_subset_varid, sample.id = shared_sampleid
  )

  colnames(genomat$genotype) <- genomat$sample.id
  rownames(genomat$genotype) <- genomat$snp.id
  cat("-- Finish convertion. \n\n")

  #### sort according to shared.sample.id ======================================
  geno_mat <- genomat$genotype[, match(shared_sampleid, genomat$sample.id)]
  saveRDS(head(geno_mat), "input_geno_mat_header.rds")

  cvrt_mat <- cvrt_dat[match(shared_sampleid, cvrt_dat$sample.id), ]
  cvrt_mat$sample.id <- NULL
  rownames(cvrt_mat) <- shared_sampleid
  saveRDS(cvrt_mat, "input_cvrt_pc.rds")
  cat("-- Finish reading covariates. \n\n")

  #### obtain phenotype data ===================================================

  cat("-- Start analyzing", chunk_name, "...\n")

  pheno_mat <- pheno_dat[match(shared_sampleid, pheno_dat$sample.id), ]
  rownames(pheno_mat) <- shared_sampleid
  pheno_mat$sample.id <- NULL
  saveRDS(pheno_mat, paste0("input_phenotype_", chunk_name, ".rds"))

  cat("-- Start matrixeQTL engine...\n")

  #### matrixeQTL ####
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

  cat("-- Finish matrixeQTL engine. Continue to format the result...\n")

  ##### QTL result -------------------------------------------------------------
  QTL_res_df <- me$all$eqtls

  if (nrow(QTL_res_df) > 0) {
    assoc_has_result <- TRUE
    cat("-- Finish QTL analysis.\n\n")

    ## rename columns; flip the direction of the effect estimate so it will be the effect of ALT
    if (model_type == "linear") {
      QTL_res_df <- QTL_res_df %>%
        mutate(statistic = (-1) * statistic, beta = (-1) * beta, FDR = NULL) %>%
        rename(variant.id = snps, phenotype = gene, beta.ca = beta, statistic.ca = statistic)
    } else if (model_type == "category") {
      c("variant.id", "phenotype", "statistic.anova", "pvalue")
      QTL_res_df <- QTL_res_df %>%
        mutate(FDR = NULL) %>%
        rename(variant.id = snps, phenotype = gene, statistic.anova = statistic)
    } else if (model_type == "interaction") {
      QTL_res_df <- QTL_res_df %>%
        mutate(statistic = (-1) * statistic, beta = (-1) * beta, FDR = NULL) %>%
        rename(variant.id = snps, phenotype = gene, beta.inter = beta, statistic.inter = statistic)
    }

    saveRDS(QTL_res_df, paste0("QTL_", gds_name, "_", chunk_name, ".rds"))

    #### draw phenotype-genotype boxplot ---------------------------------------
    if (draw_genopheno_boxplot & min(QTL_res_df$pvalue) < boxplot_p_cutoff) {
      cat("-- Start plotting genotype phenotype box plot for top significant SNPs... \n")

      ## first, obtain all SNPs passed p threshold
      QTL_res_sub <- QTL_res_df %>%
        select(variant.id, phenotype, pvalue) %>%
        filter(pvalue < boxplot_p_cutoff)

      ## subset GDS to only selected variants
      seqSetFilter(gds, variant.id = QTL_res_sub$variant.id, sample.id = shared_sampleid)

      ## obtain SNP information
      SNP_info_sub <- data.frame(
        "variant.id" = seqGetData(gds, "variant.id"),
        "REF" = seqGetData(gds, "$ref"),
        "ALT" = seqGetData(gds, "$alt"),
        "chr" = seqGetData(gds, "chromosome"),
        "pos" = seqGetData(gds, "position")
      )

      ## for each phenotype, for each chromosome, only keep the top SNP
      ## sometimes multiple SNPs have the same p-value. keep the smallest pos
      merged_sub <- merge(x = SNP_info_sub, y = QTL_res_sub, by = "variant.id", all = T)
      merged_sub_filtered <- merged_sub %>%
        group_by(phenotype, chr) %>%
        arrange(pos) %>%
        slice_min(order_by = pvalue) %>%
        distinct(phenotype, chr, .keep_all = TRUE) %>%
        as.data.frame()

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
        y = pheno_dat %>% select(sample.id, all_of(merged_sub_filtered$phenotype)),
        by = "sample.id", all.x = T, all.y = F
      )
      ## this df contains all we need but each phenotype will only be plotted against a limited set of SNPs

      ## plot for each phenotype
      for (phenotype0 in merged_sub_filtered$phenotype) {
        ## SNPs for this phenotype to plot:
        variants_to_plot <- merged_sub_filtered %>%
          filter(phenotype == phenotype0) %>%
          pull(variant.id)

        ## plot for each SNP
        for (variant0 in variants_to_plot) {
          snp_name0 <- paste0("chr", paste(merged_sub_filtered %>% filter(variant.id == variant0) %>% select(chr, pos, REF, ALT) %>% as.character(), collapse = "_"))
          snp_column0 <- as.character(variant0)
          png(paste0(phenotype0, "_", snp_name0, ".png"), width = 600, height = 400, res = 120)
          print(
            ggplot(plot_pheno_dat, aes(x = !!as.name(snp_column0), y = !!as.name(phenotype0), fill = !!as.name(snp_column0))) +
              geom_boxplot() +
              theme_minimal() +
              xlab(snp_name0) +
              theme(legend.position = "none") +
              scale_fill_manual(values = c("#E5A1C3", "#A1A1E5", "#A1E5C3")) +
              ggtitle(paste(phenotype0, "x", snp_name0))
          )
          dev.off()
        }
      }
    }

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

    cat("-- Finish saving results for", gds_name, "_", chunk_name, ".\n\n")
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
}

if (!run_assoc | !assoc_has_result) {
  if (useModel == modelLINEAR) QTL_expected_colnames <- c("variant.id", "phenotype", "statistic.ca", "pvalue", "beta.ca")
  if (useModel == modelANOVA) QTL_expected_colnames <- c("variant.id", "phenotype", "statistic.anova", "pvalue")
  if (useModel == modelLINEAR_CROSS) QTL_expected_colnames <- c("variant.id", "phenotype", "statistic.inter", "pvalue", "beta.inter")

  QTL_psedo_res <- data.frame(matrix(ncol = length(QTL_expected_colnames), nrow = 0))
  colnames(QTL_psedo_res) <- QTL_expected_colnames
  saveRDS(QTL_psedo_res, paste0("QTL_noQTL_", gds_name, "_", chunk_name, ".rds"))

  QTLcount <- data.frame(
    "phenotype" = phenonames,
    "genotype_dat" = gds_name,
    "QTL_num" = 0,
    "description" = "no_result"
  )
}

write.table(QTLcount, paste0("QTL_count_", gds_name, "_", chunk_name, ".txt"),
  quote = F, sep = "\t", col.names = T, row.names = F
)

cat("-- Finish. \n")

closefn.gds(gds)
date()
sink()

#!/usr/bin/env Rscript

## input description:
## 1. individual_phenodat_chunk  rds for each individual phenotype data
## 2. directory which has all QTL result .rds
## 3. output_result_csv  output QTL results in txt format or not
## 4. plot_mac mac threshold for plotting 
## 5. plot_resolution 
## 6. plot_size 

args <- commandArgs(trailingOnly = TRUE)
chunk_name <- gsub(".*phenodat_(.+).rds.*", "\\1", args[1])
pheno_names <- colnames(readRDS(args[1]))
pheno_names <- pheno_names[which(pheno_names != "sample.id")]

QTLres_dir <- args[2]
output_result_csv <- ifelse(args[3] == "true", TRUE, FALSE)
plot_mac <- as.numeric(args[4])
plot_resolution <- as.numeric(args[5])
plot_size <- as.numeric(args[6])

sink(paste0("5_wrap_", chunk_name, ".log"), append = FALSE, split = TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))

phenodat <- readRDS(args[1])
QTLres_files <- list.files(QTLres_dir, pattern = paste0("_", chunk_name, ".rds"))
cat("#### Fish \n")
cat(QTLres_files)
QTLres_all <- NULL

if (length(QTLres_files) == 0) {
  stop("No QTL result found in directory:", QTLres_dir)
} else if (length(QTLres_files) == 1) {
  QTLres_all <- readRDS(paste0(QTLres_dir, QTLres_files))
  cat("## Only one (1) result file found in directory:", QTLres_dir)
} else if (length(QTLres_files) > 1) {
  cat("## Found", length(QTLres_files), "QTL result files. \n")
  QTLres_all <- readRDS(paste0(QTLres_dir, QTLres_files[1]))
  for (idx in c(2:length(QTLres_files))) {
    QTLres0 <- readRDS(paste0(QTLres_dir, QTLres_files[idx]))
    QTLres_all <- rbind(QTLres_all, QTLres0)
  }
  cat("## Finish merging chunk. Output file as:", paste0("QTL_", chunk_name, ".rds"), "\n\n")
}

if (!is.null(QTLres_all)) {
  for (pheno_idx in seq(pheno_names)) {
    pheno_name0 <- pheno_names[pheno_idx]
    cat("## Start ", pheno_name0, "...\n")

    QTLres <- QTLres_all %>% filter(phenotype == pheno_name0)

    if (nrow(QTLres) > 0) {
      saveRDS(QTLres, paste0("QTL_", pheno_name0, ".rds"))
      if (output_result_csv) {
        write.csv(QTLres, paste0("QTL_", pheno_name0, ".csv"), row.names = FALSE)
      }
      cat("## Finish saving", pheno_name0, ". Output file:", paste0("QTL_", pheno_name0, ".rds/csv"), "\n")

      dir.create(paste0("plots_", pheno_name0))
      plot_dir <- paste0("plots_", pheno_name0, "/")

      #### phenotype plot ####
      png(filename = paste0(plot_dir, "hist_", pheno_name0, ".png"), width = plot_size * 1.5, height = plot_size, res = plot_resolution)
      hist(phenodat[, pheno_name0], col = "lightskyblue1", main = paste("Histogram of", pheno_name0), xlab = pheno_name0)
      dev.off()
      cat("## Finish plotting phenotype histogram. \n\n")

      #### qq plot ####
      QTLres <- QTLres %>% filter(pvalue != 0, MAC >= plot_mac)

      lambda <- median(qchisq(1 - QTLres$pvalue, 1)) / qchisq(0.5, 1)

      png(filename = paste0(plot_dir, "qq_", pheno_name0, "_lambda_", round(lambda, digits = 3), ".png"), width = plot_size * 1.5, height = plot_size * 1.5, res = plot_resolution)
      qq(QTLres$pvalue)
      dev.off()
      cat("## Finish plotting QQ plot.\n")
      cat("   Lambda:", lambda, " \n\n")

      #### manhattan plot ####
      png(filename = paste0(plot_dir, "mht_", pheno_name0, ".png"), width = plot_size * 2.5, height = plot_size * 1.5, res = plot_resolution)
      manhattan(x = QTLres, chr = "chr", bp = "pos", snp = "variant.id", p = "pvalue")
      dev.off()
      cat("## Finish plotting mht plot.\n\n")
    } else {
      cat("## No QTL found for", pheno_name0, ".\n\n")
    }
  }
}



date()
sink()

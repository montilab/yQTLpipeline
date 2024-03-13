#!/usr/bin/env Rscript

## input description:
## 1. individual_phenodat  rds for each individual phenotype data
## 2. directory which has all QTL result .rds
## 3. output_result_csv  output QTL results in txt format or not
## 4. plot_mac mac threshold for plotting
## 5. plot_resolution
## 6. plot_size

args <- commandArgs(trailingOnly = TRUE)
phenoname <- gsub(".*phenodat_(.+).rds.*", "\\1", args[1])
QTLres_dir <- args[2]
output_result_csv <- ifelse(args[3] == "true", TRUE, FALSE)
plot_mac <- as.numeric(args[4])
plot_resolution <- as.numeric(args[5])
plot_size <- as.numeric(args[6])

dir.create(paste0("plots_", phenoname))
plot_dir <- paste0("plots_", phenoname, "/")

sink(paste0(plot_dir, "5_wrap_", phenoname, ".log"), append = FALSE, split = TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(qqman))

phenodat <- readRDS(args[1])
QTLres_files <- list.files(QTLres_dir, pattern = paste0(phenoname, ".+.rds"))

QTLres <- NULL
if (length(QTLres_files) == 0) {
  stop("No QTL result found in directory:", QTLres_dir, "\n")
} else if (length(QTLres_files) == 1) {
  QTLres <- readRDS(paste0(QTLres_dir, QTLres_files))
  cat("-- Only one (1) result file found in directory:", QTLres_dir, "\n")
} else if (length(QTLres_files) > 1) {
  cat("-- Found", length(QTLres_files), "QTL result files in directory:", QTLres_dir, "\n")
  QTLres <- readRDS(paste0(QTLres_dir, QTLres_files[1]))
  for (idx in c(2:length(QTLres_files))) {
    QTLres0 <- readRDS(paste0(QTLres_dir, QTLres_files[idx]))
    QTLres <- rbind(QTLres, QTLres0)
  }
  cat("-- Finish merging.")
}

saveRDS(QTLres, paste0("QTL_", phenoname, ".rds"))
if (output_result_csv) {
  write.csv(QTLres, paste0("QTL_", phenoname, ".csv"), row.names = FALSE)
}

cat("-- Finish saving. Output file as:", paste0("QTL_", phenoname, ".rds/csv"), "\n\n")

if (!is.null(QTLres)) {
  #### Phenotype distribution plot ####
  png(filename = paste0(plot_dir, "hist_", phenoname, ".png"), width = plot_size * 1.5, height = plot_size, res = plot_resolution)
  hist(phenodat[, phenoname], col = "lightskyblue1", main = paste("Histogram of", phenoname), xlab = phenoname)
  dev.off()
  cat("-- Finish plotting phenotype histogram. \n\n")

  #### QQ plot ####
  QTLres <- QTLres %>% filter(Score.pval != 0, MAC >= plot_mac)

  lambda <- median(qchisq(1 - QTLres$Score.pval, 1)) / qchisq(0.5, 1)

  png(filename = paste0(plot_dir, "qq_", phenoname, "_lambda_", round(lambda, digits = 3), ".png"), width = plot_size * 1.5, height = plot_size * 1.5, res = plot_resolution)
  qq(QTLres$Score.pval)
  dev.off()
  cat("-- Finish plotting QQ plot.\n")
  cat("   Lambda:", lambda, " \n\n")

  #### Manhattan plot ####
  png(filename = paste0(plot_dir, "mht_", phenoname, ".png"), width = plot_size * 2.5, height = plot_size * 1.5, res = plot_resolution)
  manhattan(x = QTLres, chr = "chr", bp = "pos", snp = "variant.id", p = "Score.pval")
  dev.off()
  cat("-- Finish plotting mht plot.\n\n")

  #### Miami plot ####
  ylim <- max(-log10(QTLres$Score.pval)) + 0.2
  png(filename = paste0(plot_dir, "miami_", phenoname, ".png"), width = plot_size * 2.5, height = plot_size * 1.5, res = plot_resolution)
  par(mfrow = c(2, 1))
  par(mar = c(1.3, 3, 3, 3))
  manhattan(QTLres %>% filter(Est.ca > 0), ylim = c(0, ylim), chr = "chr", bp = "pos", snp = "variant.id", p = "Score.pval")
  par(mar = c(3, 3, 1.3, 3))
  manhattan(QTLres %>% filter(Est.ca < 0), ylim = c(ylim, 0), xlab = "", xaxt = "n", chr = "chr", bp = "pos", snp = "variant.id", p = "Score.pval")
  dev.off()
  cat("-- Finish plotting Miami plot.\n\n")
}

date()
sink()

#!/usr/bin/env Rscript

## merge SNP info obtained from GDS with QTL results ##
## input description:
##  1. gds_name  name of gds file, e.g. chromosome number
##  2. QTL_result  QTL of a phenotype of this chromosome
##  3. SNP_info  SNP information obtained from GDS

## output:
## annot_assoc_{pheno_name}_{gds_name}.rds

args <- commandArgs(trailingOnly = TRUE)

sink(paste0("4_", gsub(".*assoc_(.+).rds.*", "\\1", args[2]), "_merge_SNP_info.log"), append = FALSE, split = TRUE)
date()

gds_name <- args[1]
QTL_res <- readRDS(args[2])
SNP_info <- read.csv(args[3])
phenoname <- QTL_res$phenotype[1]

merged <- merge(x = SNP_info, y = QTL_res, by = "variant.id", all = F)

saveRDS(merged, paste0("assoc_SNPinfo_", phenoname, "_", gds_name, ".rds"))

date()
sink()

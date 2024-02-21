#!/usr/bin/env Rscript

## PCAiR ## PC_Relate ##
## input description:
##  1. gds_file  file path for .gds
##  2. snp_PCA_txtfile  "NA" if not available; simple txt file format, only one column listing the annotation id in .gds to use to calculate PC
##  3. pcair_kinthresh  threshold, 2^(-X/2), input the X here
##  4. pcair_divthresh  threshold, -2^(-X/2), input the X here
##  5. sampleid_use_file  optional user defined sample ids to use. use "NA" to include all samples in gds

## output: temp_king.rds  temp_pcair.rds  PCs.rds  temp_pcrelate.rds  GRM.rds
##         PCA_plot.png  kinship_plot.png

sink("0_PC_GRM.log", append = FALSE, split = TRUE)
date()

args <- commandArgs(trailingOnly = TRUE)
gds_file <- args[1]

if (args[2] == "NA") {
  customized_snpset <- FALSE
} else {
  customized_snpset <- TRUE
  snp_PCA_txtfile <- args[2]
}

pcair_kinthresh <- 2^(-as.numeric(args[3]) / 2)
pcair_divthresh <- -2^(-as.numeric(args[4]) / 2)

suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

#### genotype data + sample ids #### ===========================================
genodat <- seqOpen(gds_file)

genodat_sampleid <- seqGetData(genodat, "sample.id")
if (args[5] != "NA") {
  user_def_sampleid <- readLines(args[5])
  cat("Number of user defined samples:", length(user_def_sampleid), ".\n")
  sampleid_use <- intersect(user_def_sampleid, genodat_sampleid)
  cat("Number of samples available in GDS file: ", length(sampleid_use), ".\n")
} else {
  cat("Number of samples in GDS file:", length(genodat_sampleid), ".\n")
  sampleid_use <- genodat_sampleid
}

#### SNP set #### ==============================================================
if (customized_snpset) {

  ## customized SNP set ##
  ## PCAiR requires variant id in GDS file
  ## But to make the mapping easier, allow the user to customize a SNP set with annotation id
  cat("Reading SNP set...\n")
  snpset <- as.character(readLines(snp_PCA_txtfile))
  cat("Retriving SNP variant id and annotation id in merged GDS file...\n")
  SNP_IDs <- data.frame(
    "variant.id" = seqGetData(genodat, "variant.id"),
    "annot.id" = seqGetData(genodat, "annotation/id")
  )
  SNP_IDs <- SNP_IDs %>% filter(annot.id %in% snpset)
  snpset <- SNP_IDs$variant.id
  cat(
    "Finished.", length(snpset), "SNPs used.",
    length(setdiff(snpset, SNP_IDs$annot.id)), "SNPs not found in GDS. \n"
  )
} else {

  ## LD based correlation pruning ##
  cat("Start LD pruning to get SNP set for PCA...\n")
  snpset <- snpgdsLDpruning(genodat,
    method = "corr",
    slide.max.bp = 1e8, ld.threshold = sqrt(0.1)
  )
  cat("SNP set generated.\n")
}

snpset <- unlist(snpset)

#### Create SeqVarData object #### =============================================

seqData <- SeqVarData(genodat)

#### KING #### =================================================================
king <- snpgdsIBDKING(
  gdsobj = genodat, sample.id = sampleid_use,
  snp.id = snpset, verbose = TRUE
)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(king, "temp_king.rds")
cat("king saved. \n")

#### PC-AiR #### ===============================================================
pcs <- pcair(seqData,
  kinobj = kingMat, kin.thresh = pcair_kinthresh,
  divobj = kingMat, div.thresh = pcair_divthresh,
  sample.include = sampleid_use,
  snp.include = snpset,
  verbose = TRUE
)
saveRDS(pcs, "temp_pcair.rds")
pc_df <- as.data.frame(pcs$vectors)
colnames(pc_df) <- paste0("PC", 1:ncol(pcs$vectors))
pc_df$sample.id <- row.names(pcs$vectors)

#### PCA plot #### =============================================================
png("PCA_plot.png", width = 800, height = 800)
par(mfrow = c(2, 2))
plot(pc_df$PC1, pc_df$PC2, cex = 1.3, cex.axis = 1.4, cex.lab = 1.4)
plot(pc_df$PC3, pc_df$PC4, cex = 1.3, cex.axis = 1.4, cex.lab = 1.4)
plot(pc_df$PC5, pc_df$PC6, cex = 1.3, cex.axis = 1.4, cex.lab = 1.4)
plot(pc_df$PC7, pc_df$PC8, cex = 1.3, cex.axis = 1.4, cex.lab = 1.4)
dev.off()

saveRDS(pc_df, "PCs.rds")

#### PC-Relate #### ============================================================
seqSetFilter(seqData, variant.id = snpset)
iterator <- SeqVarBlockIterator(seqData, variantBlock = 20000, verbose = FALSE)
pcrel <- pcrelate(iterator,
  pcs = pcs$vectors[, 1:2],
  sample.include = sampleid_use, training.set = pcs$unrels
)

saveRDS(pcrel, "temp_pcrelate.rds")
seqResetFilter(seqData, verbose = FALSE)
cat("pc-relate saved. \n")

## kinship plot
kinship <- pcrel$kinBtwn
png("kinship_plot.png")
print(ggplot(kinship, aes(k0, kin)) +
  geom_hline(yintercept = 2^(-seq(3, 9, 2) / 2), linetype = "dashed", color = "grey") +
  geom_point(alpha = 0.5) +
  ylab("kinship estimate") +
  theme_bw())
dev.off()

## GRM
grm <- pcrelateToMatrix(pcrel, scaleKin = 2)
saveRDS(grm, "GRM.rds")
cat("GRM saved. \n")

cat("PC-AiR PC-Relate finish. \n")

date()
sink()

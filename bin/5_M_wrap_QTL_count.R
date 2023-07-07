#!/usr/bin/env Rscript

## input description:
##  1. all individual QTL count per phenotype table .txt

args <- commandArgs(trailingOnly = TRUE)

sink("4_QTL_count_wrap.log", append = FALSE, split = TRUE)
date()

cat("#### QTL per phenotype count summarize start \n")
counttable_res <- data.frame("phenotype" = character(), "QTL_num" = numeric())
for (i in seq(length(args))) {
  temp <- read.table(file = args[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counttable_res <- rbind(counttable_res, temp)
}
counttable_res <- counttable_res[complete.cases(counttable_res), ]
counttable_res <- counttable_res[order(counttable_res$QTL_num, decreasing = TRUE), ]
write.table(counttable_res, "QTL_counts_summary.txt",
  col.names = TRUE, row.names = FALSE, quote = F, sep = "\t"
)
cat("Finished. \n")


date()
sink()

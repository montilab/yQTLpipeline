dir <- gsub("example", "", getwd())

gds_list <- data.frame(
  "name" = paste0("chr", 1:22),
  "path" = paste0(dir, "data/gds/chr_", 1:22, ".gds")
)
write.table(gds_list, "../data/gds_list.csv", sep = ",", row.names = F, col.names = F)

writeLines(paste0('includeConfig \'', dir, '/configs/local.config\'

params.outdir = "$PWD/ExampleResults_G"
params.datdir = "', dir, '/data"
params.pipeline_engine = "genesis"

// Genotype data input:
params.genodat_format = "gds"
params.gds_list = "${params.datdir}/gds_list.csv"
params.genotype_or_dosage = "DS"

// Phenotype input:
params.phenotype_file = "${params.datdir}/pheno_file.csv"
params.phenotype_names = "${params.datdir}/pheno_name.txt"

// Covariates:
params.covariates_file = "${params.phenotype_file}"
params.covariates = "age,sex,PC1,PC2,PC3,PC4"
params.covariates_factor = "sex"

// PCs and GRM input:
params.PC_rds = "NA"
params.GRM_rds = "${params.datdir}/grm.rds"

// Specific SNPs and Samples to analyze:
params.snpset_assoc_txtfile = "NA"
params.userdef_sampleid_txtfile = "NA"

// optional workflow for PCA and GRM:
params.start_PC = false
params.snpset_PCA_txtfile = "NA"
params.pcair_kinthresh = 11
params.pcair_divthresh = 11

// Pipeline setups:
params.max_forks_parallel = 100
params.max_pheno_parallel = 100
params.pval_cutoff = 1
params.output_result_csv = "true"

// Plotting parameters:
params.draw_genopheno_boxplot = "true"
params.boxplot_p_cutoff = 5e-8
params.plot_mac = 3 
params.plot_resolution = 100 
params.plot_size = 400 
'),
  con = "ExampleConfig_G.config"
)

writeLines('
../nextflow -c ExampleConfig_G.config run ../Analysis.nf
../nextflow -c ExampleConfig_G.config run ../Report.nf
',
  con = "Example_G.sh"
)

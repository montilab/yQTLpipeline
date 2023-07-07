# *yQTL Pipeline* Documentations 

## 1\. Setup 

### 1\.1 Install R packages
Install the required R packages: \
``` R
install.packages("dplyr")
install.packages("tidyverse")
BiocManager::install("Biobase")
install.packages("data.table") 
BiocManager::install("SNPRelate")
BiocManager::install("SeqArray")
BiocManager::install("SeqVarTools")
install.packages("MatrixEQTL")
BiocManager::install("GENESIS")
install.packages("qqman")
```

### 1\.2 Decide which Nextflow (nf) scripts to use 
The whole pipeline is divided into three nf scripts: `Prepare.nf`, `Analysis.nf` and `Report.nf`. \
1. When the genotype data is in GDS format and genetic PCs (as well as genetic relationship matrix (GRM), when individual relatedness is present) is available, only `Analysis.nf` and `Report.nf` are needed. \
2. When the genotype data is in VCF format, and/or genetic PCs and GRM are needed, `Prepare.nf` must be run before launching `Analysis.nf` and `Report.nf`. \

The processes included in each one of the nf script: \
**Prepare.nf** \
1. `vcf_to_gds`  Convert VCF files to GDS format. \
2. `write_gds_list`  Write the resulting GDS file paths into a CSV file for Analysis.nf to use. \
3. `merge_gds`  Merge the GDS files for the following process `PCA` or `PCA_GRM`. \
4. `PCA` or `PCA_GRM`  Perform genetic PCA analysis using PC-AiR and obtain GRM using PC-Relate. \

**Analysis.nf** \
1. `get_shared_sampleid`  Obtain the shared sample ids from the input genotype files and phenotype file, and the following optional input files when provided: covariate file, PC RDS file, and user defined sample id TXT file. Those are the sample ids to be used in the analysis. \
2. `split_pheno`  Split phenotype data into chunks or single phenotype in order to analyze them in parallel. \
3. `prepare_cvrt`  Read in covariates file and process categorical covariates. \
4. `SNP_info`  Obtain SNP information from the GDS file. \
	When running GENESIS version, this includes variant.id, REF (reference allele), ALT.ca (alternative allele, coded allele), and snpID (annotation in genotype file). \
	When running MatrixeQTL version, besides the entries mentioned, chr (chromosome), pos (position), freq.ca (frequency of coded allele), freq.minor (frequency of minor allele), miss (missing rate), n.obs (number of observance), MAC (minor allele count) will also be included. GENESIS will return these information along with the association results. \
5. `QTL_analysis`  Association test of the SNPs from each of the genotype data file with each of the phenotype. \
6. `merge_info_QTL`  Combine the QTL results obtain from `QTL_analysis` with corresponding SNP information obtained from `SNP_info`. \
7. `QTL_count_wrap`  Count the number of QTLs returned by `QTL_analysis`. \

**Report.nf** \
1. `QTL_results_wrap`  For each of the phenotype, merge the QTL results from `QTL_count_wrap` in `Analysis.nf`, and draw phenotype histogram plot, QQ plot and Manhattan plot. \

### 1\.3 Download Nextflow executable
Nextflow requires a POSIX compatible system (Linux, OS X, etc.) and Java 8 (recommended Java 11 or later, up to 18) to be installed. Use the following command to download the executable. Once downloaded, optionally, you may make the nextflow file accessible by your $PATH variable so you do not need to specify the full path to nextflow each time. \

``` bash
$ curl -s https://get.nextflow.io | bash
```

Sometimes Nextflow may be already installed as a module on your machine. Use `module avail nextflow` to check. If so, you do not need to download the Nextflow executable. 

## 2\. Run pipeline
Specify input file paths and parameters in `Config.config`. See the following section for description for each input parameter. \

When running on local: \
``` bash
$ module load R
$ ./nextflow -c Config.config run Prepare.nf
$ ./nextflow -c Config.config run Analysis.nf
$ ./nextflow -c Config.config run Report.nf
```

When running on a shared computer cluster for high-performance: \
1. Set up the bash job parameters in `configs/sge.config`. \
2. Modify the 1st line of `Config.config` to use sge.config: `includeConfig 'configs/sge.config'`. \
3. Submit the code above as a bash job. \

Alternatively, when Nextflow module is already installed, use the following: \

``` bash
$ module load R
$ module load nextflow
$ nextflow -c Config.config run Prepare.nf
$ nextflow -c Config.config run Analysis.nf
$ nextflow -c Config.config run Report.nf
```


## 3\. Input files and parameters 

### 3\.1 Mandatory Inputs
- `params.outdir`  Directory path to save the results. \
- `params.pipeline_engine`  When individual relatedness is present and a genetic relationship matrix (GRM) is needed, use "genesis" or "g" to run GENESIS. Linear mixed-effect models will be used. When individual relatedness is not present and all samples are independent, use either "matrixeqtl" or "m" to run MatrixeQTL. Linear models will be used. Once this parameter is setup, all processes will start with either "G_" or "M_" indicating which engine has been selected. \
- `params.genodat_format`  Either "vcf" or "gds". \
- `params.vcf_list`  When `params.genodat_format` is "vcf", input the path to the CSV file that points to the location of the VCF files. The CSV file contains two columns without header. The first column is a user-defined name for each of the genotype file input, for example, "chr1chunk1" or "vcf1". The second column is the file path to the VCF file. \
- `params.gds_list`  When `params.genodat_format` is "gds", input the path to the CSV file that points to the location of the GDS files. The setup of the CSV file is the same as described above in `params.vcf_list`. When `params.genodat_format` is "vcf", VCF files will be converted to GDS by the pipeline using `Prepare.nf`. In this case, input "${params.outdir}/1_data/gds_list.csv" here. \
- `params.gds_list` and `params.gds_list` 
- `params.genotype_or_dosage`  Use "genotype" or "GT" to use the "GT" entry in the GDS or VCF files when doing association analysis. Alternatively, use "dosage" or "DS" to use the "DS" entry. \
- `params.phenotype_file`  Phenotype data file path. It should be a data frame in either TXT (tab seperated), CSV (comma seperated) or RDS format, with rows to be samples and columns to be phenotypes to discover QTL from, such as gene expression. It should contain an additional column "sample.id" matches the sample id in the genotype data. Optionally, it can also contain the covariates to be used. \
- `params.pheno_name` A TXT file that each row is a phenotype name to analyze. Or, input character string "all" in this parameter to treat all columns except "sample.id" in `params.phenotype_file` as phenotypes to be analyzed. \

### 3\.2 Optional Inputs
Please do *not* leave any of these optional parameters empty. Input "NA" when they are not applicable or not needed. \
- `params.covariates_file`  File path including the covariates. It should be a data frame in either TXT (tab seperated), CSV (comma seperated) or RDS format, with rows to be samples and columns to be covariates to use, such as "age", "gender" and "PC1". It should contain an additional column "sample.id" matches the sample id in the genotype data. Alternatively, this entry can be "${params.phenotype_file}" if the covariates are included in the phenotype data file. \
- `params.covariates`  A string input to specify covariates to be included. Multiple inputs are separate by comma (,) without space. Remember to also include the PCs to be used, regardless of they are provided in the `covariates_file` or they will be calculated by the pipeline. An input example would be "age,gender,facility,PC1,PC2". \
params.covariates_factor = A string input to specify covariates that are categorical variables. Multiple inputs are separate by comma (,) without space. An input example would be "gender,facility". \
- `params.PC_rds`  A data frame in RDS format that rows are samples and columns are PC names. It should contain an additional column "sample.id" to match with the sample ids in the genotype data. Alternatively, the PCs can be provided in the covariates file and this parameter can be then set to "NA". If the PCs will be provided by the pipeline using `Prepare.nf`, input "${params.outdir}/1_data/PCs.rds" here.  \
- `params.GRM_rds`  The genetic relationship matrix (GRM) in RDS format. The row and column names must match the names in genotype data. If the GRM will be provided by the pipeline using `Prepare.nf`, input "${params.outdir}/1_data/GRM.rds" here. \
- `params.snpset_assoc_txtfile`  A TXT file that each row is the ID of a SNP to be included in the analysis. The SNP IDs match the "annotation/id" entry in the genotype file. Alternatively, input "NA" here to include all SNPs available. \
- `params.userdef_sampleid_txtfile`  A TXT file for user defined subset of sample ids to use. Each row is a sample id. Alternatively, input "NA" here to include all the samples that are available. \

### 3\.3 Parameters to control analysis workflow
- `params.max_forks_parallel`  Nextflow parameter, numeric value. Specify how many processes will be running in parallel in QTL discovery step. QTL analysis may have high demend on memory and computational power, so please set this up accordingly. \
- `params.max_pheno_parallel`  Numeric value. Controls the maximum number of phenotypes to include in each of the QTL discovery process when using `MatrixeQTL` engine. `MatrixeQTL` is designed to be able to handle the association test of multiple phenotypes simultaneously, thus reduces computational time but also increases memory usage. This parameter sopecifies how many phenotypes will be input into MatrixeQTL engine at once. \
- `params.pval_cutoff`  Only QTL results below this will be saved and reported. \
- `params.output_result_csv`  Boolean, `true` or `false` (lower case, *with* quote). QTL results may be large and the pipeline saves all results in RDS format to reduce space. Set this parameter to be `true` to output the QTL results as CSV format besides RDS. \

**PCA parameters** \
- `params.start_PC`  Boolean, `true` or `false` (lower case, *without* quote). Set it to be `true` to ask the pipeline to estimate genetic PCs and GRM. \
- `params.snpset_PCA_txtfile`  When `params.start_PC = true`, the user can specify a TXT file which contains the SNP IDs they want to include to perform the PCA analysis. The SNP IDs should match the variant ID entry in the genotype file. \
- `params.pcair_kinthresh`  Numeric. The threshold value on a `kinobj`, which is the matrix of pairwise kinship coefficients, used for declaring each pair of individuals as related or unrelated. The most commonly used value is 11, declaring the threshold to be 2^(-11/2) ~ 0.022, corresponding to 4th degree relatives. See documentation for [PC-Air](https://rdrr.io/bioc/GENESIS/man/pcair.html) for more details. \
- `params.pcair_divthresh`  Numeric. Threshold value on a `divobj`, which the matrix of pairwise ancestry divergence measures, used for deciding if each pair of individuals is ancestrally divergent. The most commonly used value is 11, declaring the threshold to be -2^(-11/2) ~ -0.022. See documentation for [PC-Air](https://rdrr.io/bioc/GENESIS/man/pcair.html) for more details. \

**Plotting parameters** \
Those parameters only control the phenotype expression histograms, QQ plots, and Manhattan plots and will *not* filter the QTL result outputs in the RDS. \
- `params.plot_mac`  Numeric. The minimum minor allele count (MAC) for a SNP to be considered in the QQ plot and Manhattan plot to avoid artefects caused by too little observation. When analyzing rare SNPs, this is recommended to be 3 or 5. For common SNPs, this can be much higher, such as 20. \
- `params.plot_resolution`  Numeric. The resolution dpi to save the plots. Commonly used values including 72, 100, and 120. \
- `params.plot_size`  Numeric. Control the size of the plots (pixels). Commonly used values including 400 and 600. \


## 4\. Input examples
*Note: All TXT and CSV files must have a complete final line, i.e., an empty line at the end of the file when opened in a text editor. Otherwise, the computer can not process the text file.*

CSV file to specify genotype file paths (without header): \
|        |                                   |
|:-------|:----------------------------------|
| chr\_1 | /yQTL-Pipeline/data/gds/chr_1.gds |
| chr\_2 | /yQTL-Pipeline/data/gds/chr_2.gds |
| chr\_3 | /yQTL-Pipeline/data/gds/chr_3.gds |


Phenotype data frame with covariates: \

| sample.id | age | sex | gene1 | gene2 |
|:----------|:----|:----|------:|------:|
| HG00110   | 40  | F   | 7.001 | 6.540 |
| HG00116   | 50  | M   | 7.025 | 6.440 |
| HG00120   | 65  | F   | 8.035 | 6.420 |
| HG00128   | 35  | F   | 7.015 | 6.154 |


GRM: \

|         |    HG00110 |    HG00116 |    HG00120 |    HG00128 |
|:--------|-----------:|-----------:|-----------:|-----------:|
| HG00110 |  1.0332116 | -0.0179534 |  0.0070812 | -0.0114037 |
| HG00116 | -0.0179534 |  0.9901158 |  0.1161200 | -0.0369330 |
| HG00120 |  0.0070812 |  0.1161200 |  0.9772376 | -0.0595185 |
| HG00128 | -0.0114037 | -0.0369330 | -0.0595185 |  0.9500809 |


## 5\. Output description 
The results will be saved in the directory spcified in `params.outdir`, with the following sub directories. \

Preparations and intermediate results: \
- < 1_data >. Including covariates used, sample ids that are shared between all inputs and thus analyzed. If the input genotype data was in VCF format, the converted GDS files will also be saved here. \
- < 1_phenotype_data >. Including the phenotype data after split the original input into multiple chunks (when `params.pipeline_engine = "matrixeqtl"`) or single phenotype (when `params.pipeline_engine = "genesis"`). \
- < 2_SNP_info >. Information of each SNP obtained from GDS file. \
- < 3_individual_results >. QTL association results of each genotype data with each phenotype chunk (when `params.pipeline_engine = "matrixeqtl"`) or each phenotype (when `params.pipeline_engine = "genesis"`), without SNP information. \
- < 4_individual_results_SNPinfo >. Combined < 3_individual_results > with the corresponding SNP information in < 2_SNP_info >. \

Summary results and reports: \
- < 5_Results_Summary >. Contains the QTL results of each one of the phenotype in RDS format, the count table of the number of QTL identified, as well as the phenotype expression histogram, QQ plot and Manhattan plot of each phenotype. \
Usually, this is the folder the user would take a look at. If the QTL results are too large to be merged into a single file, you can obtain them and merge them in < 4_individual_results_SNPinfo >. \

Logs: \
All logs for each of the analyses steps will be saved in each of the folders. \

Cleaning to save space: \
*proceed with caution* Nextflow would save all the intermediate results in the work/ directory. When the pipeline successfully finishes, it's highly recommended to discard the work/ directory. \
*proceed with caution* Since QTL results can be very large, when the user is sure about the results of the association of all SNPs and all phenotypes are saved in < 5_Results_Summary >, results in < 2_SNP_info >, < 3_individual_results > and < 4_individual_results_SNPinfo > can be discarded, since < 5_Results_Summary > is simply the merged version of all of them. It's recommended to discard < 4_individual_results_SNPinfo > first since it is simply subsets of < 5_Results_Summary >. \

-end-

# *yQTL Pipeline* Documentations 

*A Nextflow- and R-based pipeline to organize the multiple steps of QTL discovery analysis. Built under Nextflow DSL2.*  

[![Built With](https://img.shields.io/badge/Built%20With-Nextflow%20DSL2-blue.svg)](https://www.nextflow.io/) 
![R](https://img.shields.io/badge/Framework-R%204.1.2-blue.svg) 
![Compatibility](https://img.shields.io/badge/Compatibility-Linux%20%2F%20OSX-orange.svg) 
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg) 
[![GitHub Issues](https://img.shields.io/github/issues/montilab/yQTL-Pipeline.svg)](https://github.com/montilab/yQTL-Pipeline/issues)

# Overview

We developed the _yQTL Pipeline_ – with ‘y’ indicating the dependent quantitative variable being modeled – to facilitate and automate large-scale QTL analysis. Prior to genome-wide association test, the pipeline supports the calculation or the direct input of pre-defined genome-wide principal components and genetic relationship matrix when applicable. User-specified covariates may also be supplied. Depending on the presence or absence of familial relatedness among the subjects, genome-wide association tests will be conducted using either a linear mixed-effect model or a linear model, respectively. Alternatively, the user has the option to treat the genotype as a categorical variable and conduct an ANOVA test, or assess the significance of the interaction between the SNP and a covariate. Through the adoption of the workflow management tool Nextflow, the pipeline parallelizes the analysis steps to optimize run-time and ensure reproducibility of the results. A user-friendly R Shiny App is also provided for the visualization of the results, including Manhattan plots of user-selected phenotype traits, and trait-QTL connection networks based on user-specified p-value thresholds. 
[[Preprint]](https://www.biorxiv.org/content/10.1101/2024.01.26.577518v1)

## Table of Content
- [1.Setup](https://github.com/montilab/yQTL-Pipeline#1-setup)  
- [2.Run the yQTL-Pipeline](https://github.com/montilab/yQTL-Pipeline#2-run-the-yqtl-pipeline)  
	- [Run Examples](https://github.com/montilab/yQTL-Pipeline#run-examples)  
- [3.Input Files and Parameters](https://github.com/montilab/yQTL-Pipeline#3-input-files-and-parameters)  
- [4.Input File Examples](https://github.com/montilab/yQTL-Pipeline#4-input-file-examples)  
- [5.Output Descriptions](https://github.com/montilab/yQTL-Pipeline#5-output-descriptions)  
- [6.Workflow Process Details](https://github.com/montilab/yQTL-Pipeline#6-workflow-process-details)  
- **[7.Shiny App](https://github.com/montilab/yQTL-Pipeline#7-Shiny-App-Usage)**  


## 1\. Setup 

### 1\.1 Install R packages
Install the required R packages:  
``` R
install.packages("dplyr")
install.packages("tidyverse")
install.packages("data.table")
install.packages("qqman")
install.packages("MatrixEQTL")
BiocManager::install("Biobase")
BiocManager::install("SNPRelate")
BiocManager::install("SeqArray")
BiocManager::install("SeqVarTools")
BiocManager::install("GENESIS")
```

### 1\.2 Decide which Nextflow (nf) scripts to use 
The whole pipeline is divided into three nf scripts: `Prepare.nf`, `Analysis.nf` and `Report.nf`.  
1. When the genotype data is in GDS format and genetic PCs (and genetic relationship matrix (GRM), when individual relatedness is present) is available, only `Analysis.nf` and `Report.nf` are needed.  
2. When the genotype data is in VCF format, and/or genetic PCs and GRM are needed, must be run `Prepare.nf` before launching `Analysis.nf` and `Report.nf`.  

See below in [6.Workflow Process Details](https://github.com/montilab/yQTL-Pipeline#6-workflow-process-details) for detailed descriptions of the processes included in each of the nf scripts.  

### 1\.3 Download the Nextflow executable
Nextflow requires a POSIX compatible system (Linux, OS X, etc.) and Java 8 (recommended Java 11 or later, up to 18) to be installed. Use the following command to download the executable. Once downloaded, optionally, you may make the nextflow file accessible by your $PATH variable so you do not need to specify the full path to nextflow each time.  

``` bash
$ curl -s https://get.nextflow.io | bash
```

Sometimes Nextflow may be already installed as a module on your machine. Use `module avail nextflow` to check. If so, you do not need to download the Nextflow executable.  

## 2\. Run the *yQTL Pipeline*  
Specify input file paths and parameters in `Config.config`. See the following [3.Input Files and Parameters](https://github.com/montilab/yQTL-Pipeline#3-input-files-and-parameters) for details.  

When running on local:  
``` bash
$ module load R/4.1.2 
$ ./nextflow -c Config.config run Prepare.nf 
$ ./nextflow -c Config.config run Analysis.nf 
$ ./nextflow -c Config.config run Report.nf 
```

When running on a shared computer cluster:  
1. Set up the bash job parameters in `configs/sge.config`.  
2. Modify the 1st line of `Config.config` to use sge.config: `includeConfig 'configs/sge.config'`.  
3. Submit the code above as a bash job.  

Alternatively, when Nextflow module is already installed, use the following:  

``` bash
$ module load R/4.1.2 
$ module load nextflow 
$ nextflow -c Config.config run Prepare.nf 
$ nextflow -c Config.config run Analysis.nf 
$ nextflow -c Config.config run Report.nf 
```

### Run Examples  
Two examples are provided in the folder < example/ >. Example input data are provided in < data/ >.

#### Example 1  
Start from VCF file input, calculate PCs, and use MatrixeQTL to obtain QTL results.  
You may also change "params.pipeline_engine" to "genesis" in the config to use GENESIS and the same example inputs in this example.  

1. In your bash, go to the directory < example/ > .  
`cd Directory/Where/You/Installed/the/Pipeline/example`  
2. Launch R and run `Rscript GenerateExample_M.R`.  
This will write out the files `ExampleConfig_M.config` and `Example_M.sh`, as well as update the < data/ > directory to setup the input file paths.  
3. Run `bash Example_M.sh`.  
4. Example results will be generated in folder < ExampleResults_M >.  

#### Example 2  
Start from GDS file input, using pre-defined PCs and GRM, and use GENESIS to obtain QTL results.  
You may also change "params.pipeline_engine" to "matrixeqtl" in the config to use MatrixeQTL and the same example inputs in this example.  

1. In your bash, go to the directory < example/ > .  
`cd Directory/Where/You/Installed/the/Pipeline/example`  
2. Launch R and run `Rscript GenerateExample_G.R`.  
This will write out the files `ExampleConfig_G.config` and `Example_G.sh`, as well as update the < data/ > directory to setup the input file paths.  
3. Run `bash Example_G.sh`.  
4. Example results will be generated in folder < ExampleResults_G >.  


## 3\. Input Files and Parameters 

### 3\.1 Mandatory inputs
- `params.outdir`  Directory to save the results.  
- `params.pipeline_engine`  When individual relatedness is present and a genetic relationship matrix (GRM) is needed, use "genesis" or "g" to run GENESIS. Linear mixed-effect models will be used. When individual relatedness is not present and all samples are independent, use "matrixeqtl" or "m" to run MatrixeQTL. All processes will start with either "G_" or "M_" indicating which engine has been selected.  
- `params.model_type`  Mandatory when running MatrixeQTL. Specify as one of the following: "linear" to run an additive linear model, "category" to treat genotype as categories and apply an ANOVA test, or "interaction" to test the interaction between SNP and the covariate specified in `params.interaction_cvrt`.  When running MatrixeQTL with model_type = "interaction", specify a covariate to test its interaction with the SNPs.  
- `params.genodat_format`  Either "vcf" or "gds".  
- `params.vcf_list` and `params.gds_list`  The path to the CSV file that points to the location of the VCF or GDS files, which contains two columns without header. The first column is the user-defined names for each of the genotype files, and the second column is the file paths. When the inputs are VCF files, must also setup `params.gds_list = "${params.outdir}/1_data/gds_list.csv"`.  
- `params.genotype_or_dosage`  Use "genotype" or "GT" to use the "GT" entry in the VCF or GDS files. Alternatively, use "dosage" or "DS" to use the "DS" entry.  
- `params.phenotype_file`  Phenotype data file path. Phenotype file is a data frame in either TXT (tab seperated), CSV (comma seperated) or RDS format, with rows to be samples and columns to be the phenotypes to discover QTL from. A column "sample.id" is required to match the samples with the genotype data. Optionally, the phenotype file can also contain the covariates to be used.  
- `params.pheno_name` A TXT file that each row is a phenotype name to be analyzed. Or, input string "all" to analyze all columns except "sample.id" in the phenotype file.  

### 3\.2 Optional inputs
Please do *not* leave any of these optional parameters empty. Input "NA" when they are not applicable or not required.  
- `params.covariates_file`  Covariates file, following the same format as the phenotype file. Alternatively, this entry can be "${params.phenotype_file}" if the covariates are included in the phenotype data file.  
- `params.covariates`  A string to specify the covariates to be included, separate by comma (,) without space. PCs to be used must be specified, regardless of whether they are already provided or they will be calculated by the pipeline. An example: "age,gender,facility,PC1,PC2".  
- `params.covariates_factor`  A string to specify which covariates are categorical variables. An example: "gender,facility".  
- `params.PC_rds`  A data frame in RDS format which rows are samples and columns are PCs. It requires a column "sample.id" which matches the sample ids in the genotype data. Alternatively, the PCs can be provided in the covariates file. If the PCs will be estimated by the pipeline, setup `params.PC_rds = "${params.outdir}/1_data/PCs.rds"`.  
- `params.GRM_rds`  The genetic relationship matrix (GRM) in RDS format, which row and column names are the sample ids in the genotype data. If the GRM will be estimated by the pipeline, setup `params.GRM_rds = "${params.outdir}/1_data/GRM.rds"`.  
- `params.snpset_assoc_txtfile`  A TXT file that each row is a SNP to be analyzed. Those SNP IDs match the "annotation/id" entry in the genotype file. Alternatively, input "NA" to include all SNPs available.  
- `params.userdef_sampleid_txtfile`  A TXT file that each row is a sample id to be included in the analysis. Alternatively, input "NA" to include all the available samples.  

### 3\.3 Parameters to control analysis workflow
- `params.max_forks_parallel`  Numeric value. Controls how many processes will be running in parallel in the QTL discovery step.   
- `params.max_pheno_parallel`  Numeric value. Controls the maximum number of phenotypes to be included in each of the QTL discovery process when using `MatrixeQTL` engine. `MatrixeQTL` is designed to handle the association test of multiple phenotypes simultaneously. Increase this parameter will result in less computational time but more memory usage.  
- `params.pval_cutoff`  Numeric value. Only QTL results below this p-value threshold will be saved and reported.  
- `params.output_result_csv`  Boolean, "true" or "false" (lower case *with* quote). By default, QTL results are saved in RDS format to save space. Set this parameter to be "true" to write out the QTL results as CSV format besides RDS.  

**PCA Parameters**  
- `params.start_PC`  Boolean, `true` or `false` (lower case, *without* quote). Set it to be `true` to estimate genetic PCs and GRM using the *yQTL Pipeline*.  
When `params.start_PC = true`, setup the following:  
- `params.snpset_PCA_txtfile`  Optionally, specify a TXT file which contains the SNP IDs to be included to obtain the genetic PCs. The SNP IDs match the annotation entry in the genotype file. Alternatively, input "NA" to prune the SNPs using LD.  
- `params.pcair_kinthresh`  Numeric value. The threshold value on a `kinobj`, which is the matrix of pairwise kinship coefficients, used for declaring each pair of individuals as related or unrelated. The most commonly used value is 11, declaring the threshold to be 2^(-11/2) ~ 0.022, corresponding to 4th degree relatives. See documentation for [PC-Air](https://rdrr.io/bioc/GENESIS/man/pcair.html) for more details.  
- `params.pcair_divthresh`  Numeric. Threshold value on a `divobj`, which the matrix of pairwise ancestry divergence measures, used for deciding if each pair of individuals is ancestrally divergent. The most commonly used value is 11, declaring the threshold to be -2^(-11/2) ~ -0.022. See documentation for [PC-Air](https://rdrr.io/bioc/GENESIS/man/pcair.html) for more details.  

**Plotting Parameters**  
These parameters only control the plots and will not filter the QTL results in RDS or CSV.  
- `params.draw_genopheno_boxplot`  Input "true" to generate genotype-phenotype boxplots of the most significant variant in each genotype file, if it passes the p-value threshold below.  
- `params.boxplot_p_cutoff`  The p-value threshold to determine if the top SNP is significant.  
- `params.plot_mac`  Numeric value. The minimum minor allele count (MAC) for a SNP to be included in the QQ plot, Manhattan plot and Miami plot to avoid artefects caused by too little observation. Commonly used values for rare SNPs is 3 or 5, while 20 for common SNPs.  
- `params.plot_resolution`  Numeric value. The resolution dpi to save the plots. Commonly used values including 72, 100, and 120.  
- `params.plot_size`  Numeric value. Control the size of the plots in pixels. Commonly used values including 400 and 600.  

**PreQC**
We have a supporting script PreQC.nf that can perform quality filter on the input VCF.  
- `params.hwe_p`  The Hardy-Weinberg Equilibrium p-value. Variants with a p-value below this threshold will be considered to violate HWE and will be removed. Input 1 to avoid filtering out any variants.  
- `params.min_mac`  Numeric. Only variants with minor allele count equal or higher than this threshold will be kept.  
- `params.max_missing_allowed`  Numeric, ranges from 0 to 1. The maximum missing rate acceptable. For example, input 0.1 will remove variants with missingness greater than 10%.  


## 4\. Input File Examples 
*Note: All TXT and CSV files must have a complete final line. i.e. You can see an empty line at the end of the file when opened in a text editor.*

CSV file to specify genotype file paths (without header):  
|        |                                   |
|:-------|:----------------------------------|
| chr\_1 | /yQTL-Pipeline/data/gds/chr_1.gds |
| chr\_2 | /yQTL-Pipeline/data/gds/chr_2.gds |
| chr\_3 | /yQTL-Pipeline/data/gds/chr_3.gds |


Phenotype data frame with covariates:  

| sample.id | age | sex | gene1 | gene2 |
|:----------|:----|:----|------:|------:|
| HG00110   | 40  | F   | 7.001 | 6.540 |
| HG00116   | 50  | M   | 7.025 | 6.440 |
| HG00120   | 65  | F   | 8.035 | 6.420 |
| HG00128   | 35  | F   | 7.015 | 6.154 |


GRM:  

|         |    HG00110 |    HG00116 |    HG00120 |    HG00128 |
|:--------|-----------:|-----------:|-----------:|-----------:|
| HG00110 |  1.0332116 | -0.0179534 |  0.0070812 | -0.0114037 |
| HG00116 | -0.0179534 |  0.9901158 |  0.1161200 | -0.0369330 |
| HG00120 |  0.0070812 |  0.1161200 |  0.9772376 | -0.0595185 |
| HG00128 | -0.0114037 | -0.0369330 | -0.0595185 |  0.9500809 |


## 5\. Output Descriptions 
The results will be saved in the directory spcified in `params.outdir`, with the following sub directories.  

### 5\.1 Preparations and intermediate results  
- < 1_data >. Including covariates and sample ids used. If the input genotype data was in VCF format, the converted GDS files will also be saved here.  
- < 1_phenotype_data >. Including the phenotype data after splitting the original input into multiple chunks (when `params.pipeline_engine = "matrixeqtl"`) or single phenotype (when `params.pipeline_engine = "genesis"`).  
- < 2_SNP_info >. Information of each SNP obtained from GDS file.  
- < 3_individual_results >. QTL association results of each genotype data with each phenotype chunk (when `params.pipeline_engine = "matrixeqtl"`) or each phenotype (when `params.pipeline_engine = "genesis"`), without SNP information.  
- < 4_individual_results_SNPinfo >. Combined < 3_individual_results > with the corresponding SNP information in < 2_SNP_info >.  

### 5\.2 Results and reports  
- **< 5_Results_Summary >**. Contains the QTL results of each one of the phenotype in RDS format, the count table of the number of QTL identified, as well as the phenotype expression histogram, QQ plot, Manhattan plot and Miami plot of each of the phenotypes. If the QTL results are too large to be merged into a single file, you can obtain them and merge them in < 4_individual_results_SNPinfo >.  

### 5\.3 Logs  
All logs for each of the analyses steps will be saved in each of the subfolders.  

### 5\.4 Cleaning to save space  
*proceed with caution*: Nextflow would save all the intermediate results in the work/ directory. When the pipeline successfully finishes, it's recommended to discard the `work/` directory.  
*proceed with caution*: QTL results can be very large. When the user is sure about the results of the association of all SNPs and all phenotypes are saved in < 5_Results_Summary >, results in < 2_SNP_info >, < 3_individual_results > and < 4_individual_results_SNPinfo > can be discarded, since < 5_Results_Summary > is the merged version of them. However, it's highly recommended to save a copy of the null-models of each of the phenotype (`3_individual_results/nullmod_[phenotype]_[gds].rds`).  

## 6\. Workflow Process Details
The processes included in each one of the nf script:  

**Prepare.nf**  
1. `vcf_to_gds`  Convert VCF files to GDS format.  
2. `write_gds_list`  Write the resulting GDS file paths into a CSV file for Analysis.nf to use.  
3. `merge_gds`  Merge the GDS files for the following process `PCA` or `PCA_GRM`.  
4. `PCA` or `PCA_GRM`  Perform genetic PCA analysis using PC-AiR and obtain GRM using PC-Relate.  

**Analysis.nf**  
1. `get_shared_sampleid`  Obtain the shared sample ids from the input genotype files and phenotype file, and the following optional input files when provided: covariate file, PC RDS file, and user defined sample id TXT file. Those are the sample ids to be used in the analysis.  
2. `split_pheno`  Split phenotype data into chunks or single phenotype in order to analyze them in parallel.  
3. `prepare_cvrt`  Read in covariates file and process categorical covariates.  
4. `SNP_info`  Obtain SNP information from the GDS file.  
	When running GENESIS version, this includes variant.id, REF (reference allele), ALT.ca (alternative allele, coded allele), and snpID (annotation in genotype file).  
	When running MatrixeQTL version, besides the entries mentioned, chr (chromosome), pos (position), freq.ca (frequency of coded allele), freq.minor (frequency of minor allele), miss (missing rate), n.obs (number of observance), MAC (minor allele count) will also be included. GENESIS will return these information along with the association results.  
5. `QTL_analysis`  Association test of the SNPs from each of the genotype data file with each of the phenotype.  
6. `merge_info_QTL`  Combine the QTL results obtain from `QTL_analysis` with corresponding SNP information obtained from `SNP_info`.  
7. `QTL_count_wrap`  Count the number of QTLs returned by `QTL_analysis`.  

**Report.nf**  
1. `QTL_results_wrap`  For each of the phenotype, merge the QTL results from `QTL_count_wrap` in `Analysis.nf`, and draw phenotype histogram plot, QQ plot, Manhattan plot and Miami plot.  


## 7\. Shiny App Usage
We developed an R Shiny App available in the App folder to facilitate downstream visualization.  
Open app.R using R studio, and click the `Run App` button at the top right.  

### 7\.1 Upload QTL analysis result as an RDS file 
![image](files:./vignette/figs/pic1_preview.png)  
The uploaded QTL result is a data frame including the following columns:  
- chr  
- pos  
- snpID  
- phenotype  
- pvalue / Score.pval  
- beta / beta.ca / Est / Est.ca
- mac / MAC  

An example is provided as `N2-acetyl,N6-methyllysine_Example.rds`.  

### 7\.2 Manhattan and Miami plots 
Select a phenotype in the dropdown menu, specify the significance p-value threshold, and provide a list of SNP to highlight in the plots (optional).  
![image](files:./vignette/figs/pic2_mht.png)  
![image](files:./vignette/figs/pic3_miami.png)  

### 7\.3 Genotype-Phenotype boxplot 
First, upload a GDS file, and a phenotype file in RDS or CSV format. Specify a variant from the "annotation/id" entry in the GDS file. Click "Retrieve Data". If successful, a new drop-down menu and plotting button will appear.  
For example, you can use 'data/gds/chr_1.gds' for the GDS file and 'data/pheno_file.csv' for the phenotype file. Select variant '1:1000156'.  
Second, select a phenotype from the new drop-down menu and click "Plot Boxplot".  
![image](files:./vignette/figs/pic6_box3.png)  

### 7\.4 Network of variants and phenotypes 
QTL results consist of multiple phenotypes and variants. Therefore, we offer an option to visualize the associations as a network. For example, upload `Network_Example.rds` to the "Choose QTL result" file input at the top.  
Next, select a p-value threshold and MAC (minor allele count) threshold, and plot the network.  

![image](files:./vignette/figs/pic7_nw1.png)  

Since genetic variants are oftentimes in high linkage disequilibrium (LD) with each other, only the top variant in each chromosome will be included.  
![image](files:./vignette/figs/pic7_nw2.png)  


-end-

#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// nextflow.enable.configProcessNamesValidation = false

log.info """\

QTL Pipeline - Analysis
==============================
Output Directory   : $params.outdir
Method             : $params.pipeline_engine
Genotype File      : $params.gds_list
Phenotype File     : $params.phenotype_file
Analyzed Phenotypes: $params.phenotype_names
Covariates         : $params.covariates
P-value Cutoff     : $params.pval_cutoff
==============================
"""

// ##### Not used process to aviod warnings from sge.config ####################
process M_vcf_to_gds {script:""" """}
process M_write_gds_list {script:""" """}
process M_merge_gds {script:""" """}
process M_PCA {script:""" """}

process M_QTL_results_wrap {script:""" """}

process G_vcf_to_gds {script:""" """}
process G_write_gds_list {script:""" """}
process G_merge_gds {script:""" """}
process G_PCA_GRM {script:""" """}

process G_QTL_results_wrap {script:""" """}

// ##### MatrixeQTL ############################################################

process M_get_shared_sampleid {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    tuple val(gds_example_name), val(gds_example)

    output:
    path '*' 
    path 'shared_sampleid.txt', emit: shared_sampleid
    
    script:
    """
    0_M_get_shared_sampleid.R $gds_example ${params.phenotype_file} ${params.covariates_file} ${params.userdef_sampleid_txtfile} 
    """
}

process M_split_pheno {
    publishDir "${params.outdir}/1_phenotype_data_chunk", mode: 'copy'
    
    input:
    path shared_sampleid
    
    output:
    path '*' 
    path 'phenodat_chunk*.rds', emit: phenodat_chunks
    
    script:
    """
    2_M_split_pheno.R ${params.phenotype_file} ${params.phenotype_names} ${params.max_pheno_parallel} $shared_sampleid 
    """
}

process M_prepare_cvrt {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    path shared_sampleid
    
    output:
    path '1_cvrt.rds', emit: cvrt_use
    path '*'
    
    script:
    """
    1_M_prepare_cvrt.R ${params.covariates_file} ${params.PC_rds} ${params.covariates} ${params.covariates_factor} $shared_sampleid 
    """
}

process M_SNP_info {
    publishDir "${params.outdir}/2_SNP_info", mode: 'copy'
    
    input:
    tuple val(gds_name), val(gds_file), file(shared_sampleid)    

    output:
    tuple val(gds_name), path('SNPinfo_*.csv'), emit: SNP_info
    path '*' 

    script:
    """
    2_M_SNP_info.R $gds_name $gds_file $shared_sampleid ${params.snpset_assoc_txtfile}
    """
}

process M_QTL_analysis {
    publishDir "${params.outdir}/3_individual_results", mode: 'copy'
    maxForks params.max_forks_parallel
    
    input:
    tuple val(gds_name), val(gds_file), file(cvrt_use), file(shared_sampleid), file(pheno_rds)

    output:
    tuple val(gds_name), path('QTL_*.rds'), emit: QTL_assoc
    path 'QTL_count_*.txt', emit: QTL_chunk_count
    path '*' 

    script:
    """
    3_M_matQTL.R $gds_name $gds_file $cvrt_use ${params.model_type} ${params.interaction_cvrt} ${params.pval_cutoff} $shared_sampleid ${params.snpset_assoc_txtfile} $pheno_rds 
    """
}

process M_merge_info_QTL {
    publishDir "${params.outdir}/4_individual_results_SNPinfo", mode: 'copy'
    
    input:
    tuple val(gds_name), file(QTL_result), file(gds_info)
    
    output:
    path 'QTL_SNPinfo_*.rds'
    path '*'

    script:
    """
    4_M_merge_info_QTL.R $gds_name $QTL_result $gds_info
    """
}

process M_QTL_count_wrap {
    publishDir "${params.outdir}/5_Results_Summary", mode: 'move'
    
    input:
    path QTL_chunk_count 
    
    output:
    path '*'

    script:
    """
    5_M_wrap_QTL_count.R $QTL_chunk_count 
    """
}


// ##### GENESIS ###############################################################

process G_get_shared_sampleid {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    tuple val(gds_example_name), val(gds_example)

    output:
    path '*' 
    path 'shared_sampleid.txt', emit: shared_sampleid
    
    script:
    """
    0_G_get_shared_sampleid.R $gds_example ${params.phenotype_file} ${params.covariates_file} ${params.userdef_sampleid_txtfile} 
    """
}

process G_split_pheno {
    publishDir "${params.outdir}/1_phenotype_data", mode: 'copy'
    
    input:
    path shared_sampleid
    
    output:
    path 'phenodat_*.rds', emit: phenodat_rds
    path '*' 
    
    script:
    """
    1_G_split_pheno.R ${params.phenotype_file} ${params.phenotype_names} $shared_sampleid 
    """
}

process G_prepare_cvrt {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    path shared_sampleid
    
    output:
    path '1_cvrt.rds', emit: cvrt_use
    path '*'
    
    script:
    """
    1_G_prepare_cvrt.R ${params.covariates_file} ${params.PC_rds} ${params.covariates} ${params.covariates_factor} $shared_sampleid 
    """
}

process G_SNP_info {
    publishDir "${params.outdir}/2_SNP_info", mode: 'copy'
    
    input:
    tuple val(gds_name), val(gds_file), file(shared_sampleid)    

    output:
    tuple val(gds_name), path('SNPinfo_*.csv'), emit: SNP_info
    path '*' 

    script:
    """
    2_G_SNP_info.R $gds_name $gds_file $shared_sampleid ${params.snpset_assoc_txtfile}
    """
}

process G_QTL_analysis {
    publishDir "${params.outdir}/3_individual_results", mode: 'copy'
    maxForks params.max_forks_parallel
    
    input:
    tuple val(gds_name), val(gds_file), file(cvrt_use), file(shared_sampleid), file(pheno_rds)

    output:
    tuple val(gds_name), path('assoc_*.rds'), emit: QTL_assoc
    path 'count_*.txt', emit: QTL_results_count
    path '*' 

    script:
    """
    3_G_nullmod_assoc.R $gds_name $gds_file ${params.genotype_or_dosage} $cvrt_use ${params.GRM_rds} $shared_sampleid ${params.snpset_assoc_txtfile} $pheno_rds ${params.pval_cutoff}
    """
}


process G_merge_info_QTL {
    publishDir "${params.outdir}/4_individual_results_SNPinfo", mode: 'copy'
    
    input:
    tuple val(gds_name), file(QTL_result), file(gds_info)
    
    output:
    path 'assoc_SNPinfo_*.rds', emit: QTL_assoc_winfo
    path '*'
    val 'fin', emit: fin

    script:
    """
    4_G_merge_info_QTL.R $gds_name $QTL_result $gds_info
    """
}

process G_QTL_count_wrap {
    publishDir "${params.outdir}/5_Results_Summary", mode: 'copy'
    
    input:
    path QTL_individual_count 
    
    output:
    path '*'

    script:
    """
    echo -e 'phenotype\tgenotype_dat\tQTL_num\tdescription' > QTL_counts_summary.txt 
    cat $QTL_individual_count >> QTL_counts_summary.txt 
    """
}


// ##### workflow ##############################################################

workflow {

genogds = Channel
    .fromPath(params.gds_list)
    .splitCsv(header: false)
    .map {row -> tuple(row[0], row[1])}

genogds0 = Channel
    .fromPath(params.gds_list)
    .splitCsv(header: false, limit: 1)
    .map {row -> tuple(row[0], row[1])}

if (params.pipeline_engine == "matrixeqtl" | params.pipeline_engine == "m") {

    snpset_assoc_txtfile = channel.from("${params.snpset_assoc_txtfile}")

    M_get_shared_sampleid(genogds0)
    shared_sampleid = M_get_shared_sampleid.out.shared_sampleid

    M_SNP_info(genogds.combine(shared_sampleid))
    
    M_split_pheno(shared_sampleid)

    M_prepare_cvrt(shared_sampleid)
    cvrt_use = M_prepare_cvrt.out.cvrt_use
    
    M_QTL_analysis(genogds.combine(cvrt_use).combine(shared_sampleid).combine(M_split_pheno.out.phenodat_chunks.flatten()))

    M_merge_info_QTL(M_QTL_analysis.out.QTL_assoc.combine(M_SNP_info.out.SNP_info, by: 0))
    
    M_QTL_count_wrap(M_QTL_analysis.out.QTL_chunk_count.collect())
    
}

if (params.pipeline_engine == "genesis" | params.pipeline_engine == "g") {
        
    G_get_shared_sampleid(genogds0)
    shared_sampleid = G_get_shared_sampleid.out.shared_sampleid
    
    G_SNP_info(genogds.combine(shared_sampleid))

    G_split_pheno(shared_sampleid)
    
    G_prepare_cvrt(shared_sampleid)
    cvrt_use = G_prepare_cvrt.out.cvrt_use    

    G_QTL_analysis(genogds.combine(cvrt_use).combine(shared_sampleid).combine(G_split_pheno.out.phenodat_rds.flatten()))

    G_merge_info_QTL(G_QTL_analysis.out.QTL_assoc.combine(G_SNP_info.out.SNP_info, by: 0))

    G_QTL_count_wrap(G_QTL_analysis.out.QTL_results_count.collect())
    
}

}


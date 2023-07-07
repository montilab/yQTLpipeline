#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// nextflow.enable.configprocessNamesValidation = false

log.info """\
-
QTL Pipeline - Report
================================
outdir    : $params.outdir
-
"""

params.output_result_csv = false

// ##### Not used process to aviod warnings from sge.config ####################
process M_vcf_to_gds {script:""" """} 
process M_write_gds_list {script:""" """} 
process M_merge_gds {script:""" """} 
process M_PCA {script:""" """} 

process M_get_shared_sampleid {script:""" """} 
process M_split_pheno {script:""" """} 
process M_prepare_cvrt {script:""" """} 
process M_SNP_info {script:""" """} 
process M_QTL_analysis {script:""" """} 
process M_merge_info_QTL {script:""" """} 
process M_QTL_count_wrap {script:""" """} 

process G_vcf_to_gds {script:""" """} 
process G_write_gds_list {script:""" """} 
process G_merge_gds {script:""" """} 
process G_PCA_GRM {script:""" """} 

process G_get_shared_sampleid {script:""" """} 
process G_split_pheno {script:""" """} 
process G_prepare_cvrt {script:""" """} 
process G_SNP_info {script:""" """} 
process G_QTL_analysis {script:""" """} 
process G_merge_info_QTL {script:""" """} 
process G_QTL_count_wrap {script:""" """} 

// ##### MatrixeQTL ############################################################

process M_QTL_results_wrap {
    publishDir "${params.outdir}/5_Results_Summary", mode: 'move'

    input:
    path phenodat_chunk
    val QTL_assoc_winfo_dir
    
    output:
    path '*'

    script:
    """
    5_M_wrap_QTL_results.R $phenodat_chunk $QTL_assoc_winfo_dir ${params.output_result_csv} ${params.plot_mac} ${params.plot_resolution} ${params.plot_size} 
    """
}


// ##### GENESIS ###############################################################

process G_QTL_results_wrap {
    publishDir "${params.outdir}/5_Results_Summary", mode: 'move'

    input:
    path individual_phenodat
    val QTL_assoc_winfo_dir
    
    output:
    path '*'

    script:
    """
    5_G_wrap_QTL_results.R $individual_phenodat $QTL_assoc_winfo_dir ${params.output_result_csv} ${params.plot_mac} ${params.plot_resolution} ${params.plot_size} 
    """
}


// ##### workflow ##############################################################

workflow { 

if (params.pipeline_engine == "matrixeqtl" | params.pipeline_engine == "m") {
    phenodat_chunk = channel.fromPath("${params.outdir}/2_phenotype_data_chunk/phenodat_*.rds").flatten()
    QTL_assoc_winfo_dir = "${params.outdir}/4_individual_results_SNPinfo/"
    M_QTL_results_wrap(phenodat_chunk, QTL_assoc_winfo_dir)
}

if (params.pipeline_engine == "genesis" | params.pipeline_engine == "g") {
    individual_phenodat = channel.fromPath("${params.outdir}/1_phenotype_data/*.rds").flatten()
    QTL_assoc_winfo_dir = "${params.outdir}/4_individual_results_SNPinfo/"
    G_QTL_results_wrap(individual_phenodat, QTL_assoc_winfo_dir)
}

}

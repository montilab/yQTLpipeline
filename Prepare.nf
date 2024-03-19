#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// nextflow.enable.configprocessNamesValidation = false

// ##### Not used process to aviod warnings from sge.config ####################

process qc_filters {script:""" """}

process M_get_shared_sampleid {script:""" """} 
process M_split_pheno {script:""" """} 
process M_prepare_cvrt {script:""" """} 
process M_SNP_info {script:""" """} 
process M_QTL_analysis {script:""" """} 
process M_merge_info_QTL {script:""" """} 
process M_QTL_count_wrap {script:""" """} 

process M_QTL_results_wrap {script:""" """}

process G_get_shared_sampleid {script:""" """} 
process G_split_pheno {script:""" """} 
process G_prepare_cvrt {script:""" """} 
process G_SNP_info {script:""" """} 
process G_QTL_analysis {script:""" """} 
process G_merge_info_QTL {script:""" """} 
process G_QTL_count_wrap {script:""" """} 

process G_QTL_results_wrap {script:""" """}

process move_boxplot {script:""" """}

// ##### MatrixeQTL ############################################################

process M_vcf_to_gds {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    tuple val(vcf_name), val(vcf_filepath)

    output:
    path '*.gds', emit: gds
    path '*'

    script:
    """
    0_M_vcf_to_gds.R $vcf_name $vcf_filepath
    """
}

process M_write_gds_list {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    val vcf_names

    output:
    path 'gds_list.csv'
    path '*'

    script:
    """
    0_M_write_gds_list.R ${params.outdir} $vcf_names 
    """    
}

process M_merge_gds {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    path gds

    output:
    path 'merged.gds', emit: merged_gds
    path '*'

    script:
    """
    0_M_merge_gds.R $gds
    """
}

process M_PCA {
    publishDir "${params.outdir}/1_data", mode: 'copy'

    input: 
    path merged_gds
    
    output:
    path '*' 

    script:
    """
    0_M_PCA.R $merged_gds ${params.snpset_PCA_txtfile} ${params.pcair_kinthresh} ${params.pcair_divthresh} ${params.userdef_sampleid_txtfile}
    """
}

// ##### GENESIS ###############################################################

process G_vcf_to_gds {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    tuple val(vcf_name), val(vcf_filepath)

    output:
    path '*.gds', emit: gds
    path '*'

    script:
    """
    0_G_vcf_to_gds.R $vcf_name $vcf_filepath
    """
}

process G_write_gds_list {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    val vcf_names

    output:
    path 'gds_list.csv'
    path '*'

    script:
    """
    0_G_write_gds_list.R ${params.outdir} $vcf_names 
    """    
}

process G_merge_gds {
    publishDir "${params.outdir}/1_data", mode: 'copy'
    
    input:
    path gds

    output:
    path 'merged.gds', emit: merged_gds
    path '*'

    script:
    """
    0_G_merge_gds.R $gds
    """
}

process G_PCA_GRM {
    publishDir "${params.outdir}/1_data", mode: 'copy'

    input: 
    path merged_gds
    
    output:
    path '*' 

    script:
    """
    0_G_PCA_GRM.R $merged_gds ${params.snpset_PCA_txtfile} ${params.pcair_kinthresh} ${params.pcair_divthresh} ${params.userdef_sampleid_txtfile}
    """
}


// ##### workflow ##############################################################

workflow { 

if (params.pipeline_engine == "matrixeqtl" | params.pipeline_engine == "m") {

    if (params.genodat_format == "vcf") {    
        genovcf = Channel
            .fromPath(params.vcf_list)
            .splitCsv(header: false)
            .map {row -> tuple(row[0], row[1])}
        M_vcf_to_gds(genovcf)
        M_write_gds_list(params.vcf_list)
    }

    if (params.start_PC) {
        if (params.genodat_format == "vcf") {
            M_merge_gds(M_vcf_to_gds.out.gds.collect())
        } 
        if (params.genodat_format == "gds") {        
            M_merge_gds(params.gds_list)
        }
        M_PCA(M_merge_gds.out.merged_gds)
    } 

}


if (params.pipeline_engine == "genesis" | params.pipeline_engine == "g") {

    if (params.genodat_format == "vcf") {    
        genovcf = Channel
            .fromPath(params.vcf_list)
            .splitCsv(header: false)
            .map {row -> tuple(row[0], row[1])}
        G_vcf_to_gds(genovcf)
        G_write_gds_list(params.vcf_list)
    }

    if (params.start_PC) {
        if (params.genodat_format == "vcf") {
            G_merge_gds(G_vcf_to_gds.out.gds.collect())
        } 
        if (params.genodat_format == "gds") {        
            G_merge_gds(params.gds_list)
        }
        G_PCA_GRM(G_merge_gds.out.merged_gds)
    } 
    
}

}

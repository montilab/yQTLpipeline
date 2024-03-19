#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

process M_QTL_results_wrap {script:""" """}

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

process G_QTL_results_wrap {script:""" """}

process move_boxplot {script:""" """}

// ##### pre-QC ############################################################

process qc_filters {
    publishDir "${params.outdir}/0_QC", mode: 'copy'
    
    input:
    tuple val(vcf_name), val(vcf_filepath)

    output:
    path '*'

    script:
    """
    vcftools \
    --gzvcf $vcf_filepath \
    --hwe ${params.hwe_p} \
    --mac ${params.min_mac} \
    --max-missing 1-${params.max_missing_allowed} \
    --recode --stdout | gzip -c > ${vcf_name}_qc1.vcf.gz

    bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' ${vcf_name}_qc1.vcf.gz -Oz -o ${vcf_name}.vcf.gz
    """
}

// ##### workflow ##############################################################

workflow { 
genovcf = Channel
    .fromPath(params.vcf_list)
    .splitCsv(header: false)
    .map {row -> tuple(row[0], row[1])}
qc_filters(genovcf)

}

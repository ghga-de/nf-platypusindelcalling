
process FILTER_BY_CRIT {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/filter' , mode: 'copy'

    input:
    tuple val(sample), file(vcfgz), file(vcf_tbi)

    output:
     tuple val(sample), path('*_postFiltered.vcf.gz'),  path('*_postFiltered.vcf.gz.tbi')   , emit: vcf
    //    path  "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_vcf =  "${sample}_postFiltered.vcf"
    def out_vcfgz =  "${sample}_postFiltered.vcf.gz"

    def filter_values = "ExAC AF ${params.crit_exac_maxmaf}"
    filter_values = filter_values + " EVS MAF ${params.crit_evs_maxmaf}"
    filter_values = filter_values + " GNOMAD_EXOMES AF ${params.crit_gnomad_exomes_maxmaf}"
    filter_values = filter_values + " GNOMAD_GENOMES AF ${params.crit_gnomad_genomes_maxmaf}"
    filter_values = filter_values + " 1K_GENOMES EUR_AF ${params.crit_1kgenomes_maxmaf}"
    filter_values = filter_values + " DBSNP CLN,COMMON nonexist,exist"
    filter_values = filter_values + " LocalControlAF_WGS AF ${params.crit_localcontrol_maxmaf}"
    filter_values = filter_values + " LocalControlAF_WES AF ${params.crit_recurrance}"
    filter_values = filter_values + " ???RECURRENCE_COL . ${params.crit_recuurenve}"
    """
    vcf_filter_by_crit.py $vcfgz $out_vcf $filter_values
    bgzip -c $out_vcf > $out_vcfgz
    tabix $out_vcfgz
    """
}

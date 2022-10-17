
process INDEL_EXTRACTION {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"

    publishDir params.outdir+'/indel_extract' , mode: 'copy'

    input:
    tuple val(meta), file(ch_vcf), file(ch_vcf_i)

    output:
    tuple val(meta), path('indel_*_somatic_functional_indels_conf_*_to_10.vcf')          , emit: somatic_functional
    tuple val(meta), path('indel_*_somatic_indels_conf_*_to_10.vcf')                     , emit: somatic_indel
    tuple val(meta), path('indel_*_somatic_ncRNA_indels_conf_*_to_10.vcf')               
    tuple val(meta), path('indel_*_germline_functional_indels_conf_*_to_10.vcf')         
    tuple val(meta), path('*.functional_var_count.txt')                                  
    path  "versions.yml"                                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def somatic_functional_indel_vcf  ="indel_${meta.id}_somatic_functional_indels_conf_${params.min_confidence_score}_to_10.vcf"

    """
    indel_extractor_v1.pl \\
        --infile=$ch_vcf \\
        --somout="indel_${meta.id}_somatic_indels_conf_${params.min_confidence_score}_to_10.vcf" \\
        --funcout=$somatic_functional_indel_vcf \\
        --ncout="indel_${meta.id}_somatic_ncRNA_indels_conf_${params.min_confidence_score}_to_10.vcf" \\
        --germlineout="indel_${meta.id}_germline_functional_indels_conf_${params.min_confidence_score}_to_10.vcf" \\
        --minconf=${params.min_confidence_score}
        
    cat $somatic_functional_indel_vcf | tail -n +2 | wc -l | cut -f1 -d " " > ${meta.id}.functional_var_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}

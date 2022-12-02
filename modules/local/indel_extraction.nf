
process INDEL_EXTRACTION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"

    input:
    tuple val(meta), file(ch_vcf), file(ch_vcf_i)

    output:
    tuple val(meta), path('indel_*_somatic_functional_indels_conf_*_to_10.vcf')          , emit: somatic_functional
    tuple val(meta), path('indel_*_somatic_indels_conf_*_to_10.vcf')                     , emit: somatic_indel
    tuple val(meta), path('indel_*_somatic_ncRNA_indels_conf_*_to_10.vcf')               
    tuple val(meta), path('indel_*_germline_functional_indels_conf_*_to_10.vcf')         
    path  "versions.yml"                                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args ?: ''

    """
    indel_extractor_v1.pl \\
        --infile=$ch_vcf \\
        --somout=indel_${prefix}_somatic_indels_conf_${params.min_confidence_score}_to_10.vcf \\
        --funcout=indel_${prefix}_somatic_functional_indels_conf_${params.min_confidence_score}_to_10.vcf \\
        --ncout=indel_${prefix}_somatic_ncRNA_indels_conf_${params.min_confidence_score}_to_10.vcf \\
        --germlineout=indel_${prefix}_germline_functional_indels_conf_${params.min_confidence_score}_to_10.vcf \\
        --minconf=${params.min_confidence_score} \\
        $args
        
    cat indel_${prefix}_somatic_functional_indels_conf_${params.min_confidence_score}_to_10.vcf | \\
    tail -n +2 | wc -l | cut -f1 -d " " > ${prefix}.functional_var_count.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}

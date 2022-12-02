// processes a report on the variants

process INDEL_JSON {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"

    
    input:
    tuple val(meta), file(vcf)

    output:
    path('*.indel.json')                       , emit: json
    path('versions.yml')                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args?: ''

    """
    indel_json_v1.0.pl $vcf > ${prefix}.indel.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS

    """
}

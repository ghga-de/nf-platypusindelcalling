// processes a report on the variants

process INDEL_JSON {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling.sif' :
    'kubran/odcf_indelcalling:v0' }"

    publishDir params.outdir+ '/screenshots' , mode: 'copy'

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)

    output:
    path('*.indel.json')                       , emit: json
    path('versions.yml')                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args?: ''

    """
    indel_json_v1.0.pl $vcf > ${meta.id}.indel.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS

    """
}

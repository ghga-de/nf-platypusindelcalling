process BCFTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_0':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_0' }"

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta), path("*.raw.vcf.gz"), path("*.raw.vcf.gz.tbi") , emit: vcf
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools \\
        sort \\
        --output indel_${prefix}.raw.vcf.gz \\
        $args \\
        $vcf

    tabix  indel_${prefix}.raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process PLATYPUS {
    tag "$meta.id"
    label 'process_intermediate'

    conda     (params.enable_conda ? "conda::platypus-variant=0.8.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/platypus-variant:0.8.1.2--py27hb698ca4_5' :
    'quay.io/biocontainers/platypus-variant:0.8.1.2--py27hb698ca4_5'}"


    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta), path('indel_*.raw.vcf.gz'), path('indel_*.raw.vcf.gz.tbi')      , emit: vcf
    tuple val(meta), path('indel_*.log')                                             , emit: platypus_log
    path  "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def bamlist = meta.iscontrol == '1' ? "${control},${tumor}" : "${tumor}"

    """
    platypus callVariants \\
        --nCPU=${task.cpus}\\
        --bamFiles=$bamlist \\
        --output=indel_${prefix}.raw.vcf \\
        --refFile=$ref \\
        --logFileName=indel_${prefix}.log \\
        $args

    bgzip  --threads ${task.cpus} -c indel_${prefix}.raw.vcf > indel_${prefix}.raw.vcf.gz
    tabix  indel_${prefix}.raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platypus: 0.8.1.1.3
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """

}

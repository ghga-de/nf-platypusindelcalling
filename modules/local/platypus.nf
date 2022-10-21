
process PLATYPUS {
    tag "$meta.id"
    label 'process_intermediate'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'platypus_0.8.1.1.3.sif' :
    'kubran/platypus:0.8.1.1-3'}"

 //   publishDir params.outdir , mode: 'copy'

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta), path('indel_*.vcf.gz'), path('indel_*.vcf.gz.tbi')      , emit: vcf
    tuple val(meta), path('indel_*.vcf')                                     
    tuple val(meta), path('indel_*.log')                                     , emit: platypus_log
    path  "versions.yml"                                                     , emit: versions

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
        --output=indel_${prefix}.vcf \\
        --refFile=$ref \\
        --logFileName=indel_${prefix}.log \\
        $args

    bgzip  --threads ${task.cpus} -c indel_${prefix}.vcf > indel_${prefix}.vcf.gz
    tabix  indel_${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platypus: 0.8.1.1.3
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """

}


process PLATYPUS {
    tag "$meta.id"
    label 'process_intermediate'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'platypus_0.8.1.1.3.sif' :
    'kubran/platypus:0.8.1.1-3'}"

    publishDir params.outdir , mode: 'copy'

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta), path('indel_*.vcf.gz'), path('indel_*.vcf.gz.tbi')      , emit: vcf
    tuple val(meta), path('indel_*.vcf')                                        , emit: vcf_temp
    tuple val(meta), path('indel_*.log')                                        , emit: platypus_log
    path  "versions.yml"                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:

    def out_vcf     = "indel_${meta.id}.vcf"
    def out_vcfgz   = "indel_${meta.id}.vcf.gz"
    def out_log     = "indel_${meta.id}.log"

    def bamlist = meta.iscontrol == '1' ? "${control},${tumor}" : "${tumor}"


    """
    platypus callVariants \\
        --nCPU=${params.max_cpus} \\
        --bamFiles=$bamlist \\
        --output=$out_vcf \\
        --refFile=$ref \\
        --logFileName=$out_log \\
        ${params.platypus_params}

    bgzip  --threads ${task.cpus} -c $out_vcf > $out_vcfgz
    tabix $out_vcfgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platypus: 0.8.1.1.3
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """

}

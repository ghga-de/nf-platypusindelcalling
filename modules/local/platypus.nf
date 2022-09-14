
process PLATYPUS {
    tag "$meta.id"
    label 'process_intermediate'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'platypus_0.8.1.1.3.sif' :
    'kubran/platypus:0.8.1.1-3'}"

    publishDir params.outdir+'/vcf' , mode: 'copy'

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta), path('*.platypus.vcf.gz'), path('*.platypus.vcf.gz.tbi')      , emit: vcf
    tuple val(meta), path('*.platypus.vcf')                                        , emit: vcf_temp
    tuple val(meta), path('*.platypus.log')                                        , emit: platypus_log
    path  "versions.yml"                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:

    def out_vcf     = "${meta.id}.platypus.vcf"
    def out_vcfgz   = "${meta.id}.platypus.vcf.gz"
    def out_log     = "${meta.id}.platypus.log"

    def bamlist = meta.iscontrol == '1' ? "${control},${tumor}" : "${tumor}"


    """
    platypus callVariants \\
        --nCPU=${params.max_cpus} \\
        --bamFiles=$bamlist \\
        --output=$out_vcf \\
        --refFile=$ref \\
        --logFileName=$out_log \\
        ${params.optimized}

    bgzip  --threads ${task.cpus} -c $out_vcf > $out_vcfgz
    tabix $out_vcfgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    platypus: 0.8.1.1.3 | sed 's/^.*platypus //; s/Using.*\$//')
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS

    """

}

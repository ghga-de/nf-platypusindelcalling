
process PLATYPUS {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/environment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'shub://IARCbioinfo/platypus-nf' :
        'iarcbioinfo/platypus-nf' }"

    publishDir '${params.outdir}' , mode: 'copy', pattern: '*.vcf'


    input:
    tuple val(sample), path(tumor), path(tumor_bai), path(control), path(control_bai)
    path(ref)
    path(ref_fai)

    output:
    tuple val(sample), path('*platypus.vcf')                      , emit: platypus_vcf
    path  "versions.yml"                                          , emit: versions

    script:

    def opt_options = "--badReadsThreshold=0 --qdThreshold=0 --rmsmqThreshold=20 --hapScoreThreshold=10 --scThreshold=0.99"
    def options_arg = "${params.optimized}" == "" ? "${params.options}" : "${opt_options}"
    def out_vcf = "${sample}.platypus.vcf"

    """

    platypus callVariants \\
    --nCPU=${params.max_cpus}  \\
    --bamFiles=${tumor} \\
    --output='file.vcf' \\
    --refFile=${ref} \\
    $options_arg

    sed 's/^##FORMAT=<ID=NV,Number=.,/##FORMAT=<ID=NV,Number=A,/1g' file.vcf | sed 's/^##FORMAT=<ID=NR,Number=.,/##FORMAT=<ID=NR,Number=A,/1g' | sed 's/^##INFO=<ID=FR,Number=.,/##INFO=<ID=FR,Number=A,/1g' > $out_vcf
    rm 'file.vcf'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    platypus: \$(echo \$(platypus --version 2>&1) | sed 's/^.*platypus //; s/Using.*\$//')
    END_VERSIONS
    """
}

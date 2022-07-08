
process PLATYPUS {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/environment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/platypus:0.8.1.1-3' }"

    publishDir params.outdir+'/platypus' , mode: 'copy'


    input:
    tuple val(meta), path(tumor), path(control)
    tuple path(ref), path(ref_fai)

    output:
    tuple val(meta), path('*.platypus.vcf')                     , emit: vcf
    tuple val(meta), path('*.platypus.log')                     , emit: platypus_log
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

script:

    def opt_options = " --genIndels=1 --genSNPs=0 --verbosity=1 --bufferSize=100000 --maxReads=5000000 --minFlank=0"
    def options_arg = "${params.optimized}" == "" ? "${params.options}" : "${opt_options}"
    def out_vcf     = "${meta.id}.platypus.vcf"
    def out_log     = "${meta.id}.platypus.log"


    if (meta.iscontrol == 1)
        {
        """
        platypus callVariants \\
        --nCPU=${params.max_cpus}  \\
        --bamFiles=$control[0],$tumor[0] \\
        --output=$out_vcf \\
        --refFile=$ref \\
        --logFileName=$out_log \\
        $options_arg

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        platypus: \$(echo \$(platypus --version 2>&1) | sed 's/^.*platypus //; s/Using.*\$//')
        END_VERSIONS

        """
        }

    else {
        """
        platypus callVariants \\
        --nCPU=${params.max_cpus}  \\
        --bamFiles=$tumor[0] \\
        --output=$out_vcf \\
        --refFile=$ref \\
        --logFileName=$out_log \\
        $options_arg

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        platypus: 0.8.1.1-3 | sed 's/^.*platypus //; s/Using.*\$//')
        END_VERSIONS
         """
        }

}

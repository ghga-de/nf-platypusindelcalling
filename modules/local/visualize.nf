
process VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/visualize' , mode: 'copy'

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple path(ref), path(ref_fai)
    tuple path(repeatmasker), path(repeatmasker_tbi)

    output:
    tuple val(meta), path('*.pdf')                   , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:

    if (meta.iscontrol == '1') {
        """
        visualize.py  \\
        --vcf=$vcfgz  \\
        --control=$meta.control_bam  \\
        --tumor=$meta.tumor_bam  \\
        --ref=$ref  \\
        --prefix=indel_  \\
        --window=${params.window_size}  \\
        --annotations=$repeatmasker
        """
        }
    else {
        """
        visualize.py  \\
        --vcf=$vcfgz  \\
        --tumor=$meta.tumor_bam  \\
        --ref=$ref  \\
        --prefix=indel_  \\
        --window=${params.window_size}  \\
        --annotations=$repeatmasker
        """
         }
}

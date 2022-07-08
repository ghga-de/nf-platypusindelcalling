
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
    file (ref)
    file (repeatmasker)

    output:
    tuple val(meta), path('*.pdf')   , emit: pdfs
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    if (meta.iscontrol == 1) {
        """
        visualize_somatic.sh -i $vcfgz -c false -r $ref -n $meta.control_bam -t $meta.tumor_bam -w ${params.window_size} -a $repeatmasker -m ${params.max_var_screenshots}

        """
        }
    else {
        """
        visualize_somatic.sh -i $vcfgz -c true -r $ref -n 'null' -t $meta.tumor_bam -w ${params.window_size} -a $repeatmasker -m ${params.max_var_screenshots}
        """
    }
}

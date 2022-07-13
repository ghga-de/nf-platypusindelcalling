
process CHECK_VARIANTS_SIZE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/perl_annotate:v0' }"
    publishDir params.outdir+'/currupted' , mode: 'copy'

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)

    output:

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args?: ''

    """
    check_variants_size.sh -i $vcf -m ${params.max_var_screenshots}
    """
}

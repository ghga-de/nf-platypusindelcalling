// samtools view will be used to get sample names
process CHECK_IF_CORRUPTED {
    tag "$sample"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/perl_annotate:v0' }"
    publishDir params.outdir+'/curropted' , mode: 'copy'
    input:
    tuple val(sample), file(file)
    val(extraVCFColumn)

    output:
    path("corrupted.txt")                  , emit: corrupted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    grep -v "^#" ${file} | cut -f ${extraVCFColumn} | sort | uniq | grep -v "^\$" | wc -l > corrupted.txt

    if [[ \$(cat corrupted.txt) -gt 0 ]]
    then
        echo "VCF is corrupted"

    """

}

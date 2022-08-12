// samtools view will be used to get sample names
process CHECK_IF_CORRUPTED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling.sif' :
    'kubran/odcf_indelcalling:v0' }"


    input:
    tuple val(meta), file(vcf), file(vcf_tbi)

    output:
 //   path("corrupted.txt")                  , emit: corrupted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args?: ''
    def nocontrol = meta.iscontrol == '1' ? 'false': 'true'
    
// If there is no control (iscontrol=1), isNoControlWorkflow arg is false

    """
    corrupted.sh -i $vcf -c $nocontrol
    """
    
}

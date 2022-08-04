// samtools view will be used to get sample names
process CHECK_IF_CORRUPTED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'runtinda_perkenv.2.siff' :
    'kubran/odcf_indelcalling:v0' }"
    publishDir params.outdir+'/currupted' , mode: 'copy'

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)

    output:
 //   path("corrupted.txt")                  , emit: corrupted

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args?: ''

// If there is no control (iscontrol=1), isNoControlWorkflow arg is false
    if (meta.iscontrol == '1'){
        """
        corrupted.sh -i $vcf -c false
        """
    }
    else {
        """
        corrupted.sh -i $vcf -c true
        """
    }
}

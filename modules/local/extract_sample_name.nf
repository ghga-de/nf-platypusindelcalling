
process EXTRACT_SAMPLE_NAME {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/confidenceannotation' , mode: 'copy'

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control), path(control_bai)

    output:
    tuple val(meta), path('*.tumorname.txt')     , emit: tumor_name
    tuple val(meta), path('*.controlname.txt')   , emit: control_name
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    if (meta.iscontrol == '1')
    {
        """
        samtools view -H $meta.control_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.controlname.txt
        samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.tumorname.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else
    {
        """
        samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.tumorname.txt
        touch ${meta.id}.controlname.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}

process CREATE_CONTIGHEADER {
    tag "$meta.id"
    label 'process_single'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v1' :'kubran/odcf_platypusindelcalling:v1' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    tuple val(meta), path ("*.header")    , emit: header
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    samtools view -H $tumor > header.txt

    awk -F '[:\\t ]+' '/^@SQ/ {print "##contig=<ID=" \$3 ",length=" \$5 ">"}' header.txt > ${meta.id}.header 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
process GREP_SAMPLENAME {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    tuple val(meta), env(tumorname), env(controlname)   , emit: samplenames
    path "versions.yml"                                 , emit: versions

    script:     
    if (meta.iscontrol == 1) {
        """
        controlname=`samtools view -H $control | grep '^@RG' | sed 's/.*SM:\\([^[:space:]]*\\).*/\\1/' | uniq`
        tumorname=`samtools view -H $tumor | grep '^@RG' | sed 's/.*SM:\\([^[:space:]]*\\).*/\\1/' | uniq`

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        controlname='dummy'
        tumorname=`samtools view -H $tumor | grep '^@RG' | sed 's/.*SM:\\([^[:space:]]*\\).*/\\1/' | uniq`
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
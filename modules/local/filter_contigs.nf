process FILTER_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    conda (params.enable_conda ? "${moduleDir}/environment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'quay.io/biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(index), path(contigs_bed)

    output:
    tuple val(meta), path("*_contigsfiltered.vcf.gz"), path("*_contigsfiltered.vcf.gz.tbi"), emit: filtered_vcf
    path "versions.yml"                                                                    , emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    
    """
    # Convert the 1-column contig list into a valid 3-column BED file
    awk '{print \$1 "\\t0\\t500000000"}' ${contigs_bed} > targets.bed

    bcftools view \\
        --targets-file ^targets.bed \\
        -O u \\
        ${vcf} | \\
    bcftools sort \\
        --temp-dir . \\
        -O z \\
        -o ${prefix}_contigsfiltered.vcf.gz

    tabix -p vcf ${prefix}_contigsfiltered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( (bcftools --version 2>&1 | head -n1) || (bcftools 2>&1 | head -n1) | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
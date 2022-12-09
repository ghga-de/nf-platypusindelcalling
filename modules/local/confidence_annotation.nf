process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple val(meta), file(a), file(b), val(tumorname), val(controlname)

    output:
    tuple val(meta), path('*.vcf.gz') ,   path('*.vcf.gz.tbi')    , emit: vcf_ann
    path  "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def samples    = meta.iscontrol == "1" ? "--controlColName=$controlname --tumorColName=$tumorname" : "--nocontrol --tumorColName=$tumorname"

    """
    confidenceAnnotation_Indels.py --infile=$vcfgz \\
        $samples $args | tee indel_${prefix}.ann.vcf | cut -f 1-11 > indel_${prefix}.conf.vcf

    bgzip -c indel_${prefix}.ann.vcf > indel_${prefix}.vcf.gz
    tabix indel_${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """
    
}

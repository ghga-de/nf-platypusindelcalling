process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v7.sif' :'kubran/odcf_indelcalling:v7' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple val(meta), file(a), file(b), val(tumorname), val(controlname)

    output:
    tuple val(meta), path('*.ann.vcf.gz') ,   path('*.ann.vcf.gz.tbi')    , emit: vcf_ann
    path  "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def samples    = meta.iscontrol == "1" ? "--controlColName=$controlname --tumorColName=$tumorname" : "--nocontrol --tumorColName=$tumorname"

    """
    confidenceAnnotation_Indels.py --infile=$vcfgz --skip_order_check \\
    $samples $args | tee ${prefix}.ann.vcf | cut -f 1-11 > ${prefix}.conf.vcf

    bgzip -c ${prefix}.ann.vcf > ${prefix}.ann.vcf.gz
    tabix ${prefix}.ann.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """
    
}

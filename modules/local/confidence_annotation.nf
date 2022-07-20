### Get the real names of the columns created by Platypus
process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/confidenceannotation' , mode: 'copy'

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple val(meta), file(tumor_name)
    tuple val(meta), file(control_name)

    output:
    tuple val(meta), path('*.conf.vcf.gz'),  path('*.conf.vcf.gz.tbi')   , emit: vcf_conf
    tuple val(meta), path('*.ann.vcf.gz'),   path('*.ann.vcf.gz.tbi')    , emit: vcf
    path  "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def temp_vcf   = "${meta.id}.ann.vcf"
    def temp_vcfgz = "${meta.id}.ann.vcf.gz"
    def out_vcf    = "${meta.id}.conf.vcf"
    def out_vcfgz  = "${meta.id}.conf.vcf.gz"

    if (meta.iscontrol == '1')
    {
        """
        confidenceAnnotation_Indels.py --infile=$vcfgz --skip_order_check --controlColName=\$(cat $control_name) --tumorColName=\$(cat $tumor_name) \\
        ${params.confidence_opts_indel} | tee $temp_vcf | cut -f 1-11 > $out_vcf

        bgzip -c $out_vcf > $out_vcfgz
        tabix $out_vcfgz

        bgzip -c $temp_vcf > $temp_vcfgz
        tabix $temp_vcfgz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
        END_VERSIONS

        """
    }
    else {
        """
        confidenceAnnotation_Indels.py --infile=$vcfgz  --skip_order_check --nocontrol --tumorColName=\$(cat $tumor_name) \\
        ${params.confidence_opts_indel} | tee $temp_vcf | cut -f 1-11 > $out_vcf

        bgzip -c $out_vcf > $out_vcfgz
        tabix $out_vcfgz

        bgzip -c $temp_vcf > $temp_vcfgz
        tabix $temp_vcfgz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
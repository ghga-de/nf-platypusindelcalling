
process CONFIDENCEANNOTATION {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/confidenceannotation' , mode: 'copy'

    input:
    tuple val(sample), file(vcfgz), file(vcf_tbi)
    tuple val(sample), file(tumor), file(tumor_bai), file(control), file(control_bai), val(iscontrol)

    output:
    tuple val(sample), path('*.conf.vcf.gz'),  path('*.conf.vcf.gz.tbi')   , emit: vcf_conf
    tuple val(sample), path('*.ann.vcf.gz'),  path('*.ann.vcf.gz.tbi')     , emit: vcf
    path  "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def temp_vcf =   "${sample}.ann.vcf"
    def temp_vcfgz =   "${sample}.ann.vcf.gz"
    def out_vcf =    "${sample}.conf.vcf"
    def out_vcfgz =  "${sample}.conf.vcf.gz"

    if (iscontrol == 1)
    {
        """
        samtools view -H $control | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > controlname.txt
        samtools view -H $tumor | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > tumorname.txt
        confidenceAnnotation_Indels.py --infile=$vcfgz --controlColName=\$(cat controlname.txt) --tumorColName=\$(cat tumorname.txt) \\
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
        samtools view -H $tumor | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > tumorname.txt
        confidenceAnnotation_Indels.py --infile=$vcfgz --nocontrol --tumorColName=\$(cat tumorname.txt) \\
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

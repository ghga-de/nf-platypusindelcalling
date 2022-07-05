
process CONFIDENCEANNOTATION_NOCONTROL {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/python2.7' }"

    publishDir params.outdir+'/confidenceannotation' , mode: 'copy'

    input:
    tuple val(sample), file(vcfgz), file(vcf_tbi)
    tuple val(sample), file(tumor), file(tumor_bai), val(iscontrol)

    output:
    tuple val(sample), path('*.conf.vcf.gz'),  path('*.conf.vcf.gz.tbi')   , emit: vcf
    //    path  "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def temp_vcf = "${sample}.annotated.vcf"
    def out_vcf =  "${sample}.conf.vcf"
    def out_vcfgz =  "${sample}.conf.vcf.gz"
    // if there is no control you will do something else here
    """
    samtools view -H $tumor | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > tumorname.txt
    confidenceAnnotation_Indels.py --infile=$vcfgz --nocontrol --tumorColName=\$(cat tumorname.txt) \\
    ${params.confidence_opts_indel} | tee $temp_vcf | cut -f 1-11 > $out_vcf

    bgzip -c $out_vcf > $out_vcfgz
    tabix $out_vcfgz
    """
}

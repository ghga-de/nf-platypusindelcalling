// vcf_pipeAnnotator.sh
// Deep annotation
process PIPE_ANNOTATOR {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling.sif' :
    'kubran/odcf_indelcalling:v0' }"

    publishDir params.outdir+ '/${meta.id}'+'/annotate_vcf' , mode: 'copy'

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)
    tuple path(enchangers), path(enchangers_i)
    tuple path(cpgislands), path(cpgislands_i)
    tuple path(tfbscons), path(tfbscons_i)
    tuple path(encode_dnase), path(encode_dnase_i)
    tuple path(mirnas_snornas), path(mirnas_snornas_i)
    tuple path(cosmic), path(cosmic_i)
    tuple path(mirbase), path(mirbase_i)
    tuple path(mir_targets), path(mir_targets_i)
    tuple path(cgi_mountains), path(cgi_mountains_i)
    tuple path(phastconselem), path(phastconselem_i)
    tuple path(encode_tfbs), path(encode_tfbs_i)

    output:
    tuple val(meta), path('*.deepanno.vcf.gz'), path('*.deepanno.vcf.gz.tbi') , emit: vcf
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    zcat < $vcf | \\
    annotate_vcf.pl -a - -b $enchangers --bFileType=bed --columnName='Enhancers' | \\
    annotate_vcf.pl -a - -b $cpgislands --bFileType=bed --columnName='CpGislands' | \\
    annotate_vcf.pl -a - -b $tfbscons --bFileType=bed --columnName='TFBScons' | \\
    annotate_vcf.pl -a - -b $encode_dnase --bFileType=bed --columnName='ENCODE_DNASE' | \\
    annotate_vcf.pl -a - -b $mirnas_snornas --bFileType=bed --columnName='miRNAs_snoRNAs' | \\
    annotate_vcf.pl -a - -b $mirbase --bFileType=bed --columnName='miRBase18' | \\
    annotate_vcf.pl -a - -b $cosmic --bFileType=bed --columnName='COSMIC' --bAdditionalColumns=7,8,9 --reportLevel=1 | \\
    annotate_vcf.pl -a - -b $mir_targets --bFileType=bed --columnName='miRNAtargets' | \\
    annotate_vcf.pl -a - -b $cgi_mountains --bFileType=bed --columnName='CgiMountains' --bAdditionalColumns=4 | \\
    annotate_vcf.pl -a - -b $phastconselem --bFileType=bed --columnName='phastConsElem20bp' --bAdditionalColumns=4 | \\
    annotate_vcf.pl -a - -b $encode_tfbs --bFileType=bed --columnName='ENCODE_TFBS' > ${meta.id}.deepanno.vcf

    bgzip -c ${meta.id}.deepanno.vcf > ${meta.id}.deepanno.vcf.gz
    tabix ${meta.id}.deepanno.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS
    """
}

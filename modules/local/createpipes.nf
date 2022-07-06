// PROCESS ANNOVAR table_annovar with UCSC files

process CREATEPIPES {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/perl_annotate:v0' }"

    //conda wont work here, database files are embeded inside the docker file
    //'noelnamai/annovar:4.18' is also a docker from https://github.com/noelnamai/aws-mutation-calling/blob/master/main.nf

    publishDir params.outdir+'/createpipes' , mode: 'copy'

    input:
    tuple val(sample), file(ch_vcf), file(ch_vcf_i)
    tuple path(repeatmasker), path(repeatmasker_i)
    tuple path(dacblacklist), path(dacblacklist_i)
    tuple path(dukeexcluded), path(dukeexcluded_i)
    tuple path(hiseqdepth), path(hiseqdepth_i)
    tuple path(selfchain), path(selfchain_i)
    tuple path(mapability), path(mapability_i)
    tuple path(simpletandemrepeats), path(simpletandemrepeats_i)

    output:
    tuple val(sample), path('*.annotated.vcf.gz'), path('*.annotated.vcf.gz.tbi')   , emit: vcf
    path  "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def tempname = "${sample}.annotated.vcf"
    def tempnamegz = "${sample}.annotated.vcf.gz"
    """
    zcat < $ch_vcf | \\
    annotate_vcf.pl -a - -b $repeatmasker --bFileType=bed --columnName='REPEAT_MASKER' | \\
    annotate_vcf.pl -a - -b $dacblacklist --bFileType=bed --columnName='DAC_BLACKLIST' | \\
    annotate_vcf.pl -a - -b $dukeexcluded --bFileType=bed --columnName='DUKE_EXCLUDED' | \\
    annotate_vcf.pl -a - -b $hiseqdepth --bFileType=bed --columnName='HISEQDEPTH' | \\
    annotate_vcf.pl -a - -b $selfchain --bFileType=bed --columnName='SELFCHAIN' --bAdditionalColumns=4 --maxNrOfMatches=5 | \\
    annotate_vcf.pl -a - -b $mapability --bFileType=bed --columnName='MAPABILITY' --breakPointMode --aEndOffset=1 | \\
    annotate_vcf.pl -a - -b $simpletandemrepeats --bFileType=bed --columnName='SIMPLE_TANDEMREPEATS' --bAdditionalColumns=4 > $tempname

    bgzip -c $tempname > $tempnamegz
    tabix $tempnamegz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
    END_VERSIONS
    """
}

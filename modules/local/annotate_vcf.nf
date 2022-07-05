
process ANNOTATE_VCF {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/perl_annotate:v0' }"

    publishDir params.outdir+'/annotate' , mode: 'copy'

    input:
    tuple val(sample), file(vcf), file(vcf_tbi)
    tuple path(kgenome), path(kgenome_i)
    tuple path(dbsnpindel), path(dbsnpindel_i)
    tuple path(dbsnpsnv), path(dbsnpsnv_i)
    tuple path(exac), path(exac_i)
    tuple path(evs), path(evs_i)
    tuple path(localcontrolwgs), path(localcontrolwgs_i)
    tuple path(localcontrolwes), path(localcontrolwes_i)
    tuple path(gnomedgenomes), path(gnomedgenomes_i)
    tuple path(gnomedexomes), path(gnomedexomes_i)

    output:
    tuple val(sample), path('*.ForAnnovar.bed')                    , emit: forannovar
    tuple val(sample), path('*.vcf')                               , emit: unziped_vcf
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def for_annovar = "${sample}.ForAnnovar.bed"
    def temp_name = "${sample}.tmp"
    def out_vcf =  "${sample}.vcf"

    """
    zcat < $vcf | \\
    annotate_vcf.pl -a - -b $dbsnpindel --columnName='DBSNP' --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=10 --minOverlapFraction=0.7 --maxBorderDistanceSum=20 --maxNrOfMatches=5 | \\
    annotate_vcf.pl -a - -b $kgenome --columnName='1K_GENOMES' --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=10 --minOverlapFraction=0.7 --maxBorderDistanceSum=20 --maxNrOfMatches=5 | \\
    annotate_vcf.pl -a - -b $exac --columnName='ExAC' --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $evs --columnName='EVS' --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $gnomedexomes --columnName='GNOMAD_EXOMES' --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $gnomedgenomes --columnName='GNOMAD_GENOMES' --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $localcontrolwgs --columnName='LocalControlAF_WGS' --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $localcontrolwes --columnName='LocalControlAF_WES' --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType | \\
    tee $temp_name | vcf_to_annovar.pl ${params.chr_prefix} ${params.chr_suffix} > $for_annovar
    mv $temp_name $out_vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS
    """
}

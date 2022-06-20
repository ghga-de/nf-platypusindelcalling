
process ANNOTATE_VCF {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    '' }"

    publishDir params.outdir+'/annotate' , mode: 'copy'

    input:
    tuple val(sample), file(vcf), file(vcf_tbi)
    path(kgenome)
    path(dbsnpindel)
    path(dbsnpsnv)
    path(exac)
    path(evs)
    path(localcontrolwgs)
    path(localcontrolwes)
    path(gnomedgenomes)
    path(gnomedexomes)

    output:
    tuple val(sample), path('*.ForAnnovar.bed')                    , emit: forannovar
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def for_annovar = "${sample}.ForAnnovar.bed"

    """
    zcat $vcf | \\
    annotate_vcf.pl -a - -b $dbsnpindel --columnName=DBSNP --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=10 --minOverlapFraction=0.7 --maxBorderDistanceSum=20 --maxNrOfMatches=5 | \\
    annotate_vcf.pl -a - -b $kgenome --columnName=1K_GENOMES --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=10 --minOverlapFraction=0.7 --maxBorderDistanceSum=20 --maxNrOfMatches=5 | \\
    annotate_vcf.pl -a - -b $exac --columnName=ExAC --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $evs --columnName=EVS --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $gnomedexomes --columnName=GNOMAD_EXOMES --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $gnomedgenomes --columnName=GNOMAD_GENOMES --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $localcontrolwgs --columnName=LocalControlAF_WGS --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType \\
    annotate_vcf.pl -a - -b $localcontrolwes --columnName=LocalControlAF_WES --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType \\
    tee temprorary | annotate_vcf.pl ${params.chr_prefix} params.chr_prefix > $for_annovar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS
    """
}

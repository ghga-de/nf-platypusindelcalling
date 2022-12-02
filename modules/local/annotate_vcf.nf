//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    ' docker://kubran/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"

    input:
    tuple val(meta)            , file(vcf)     , file(vcf_tbi), val(tumorname), val(controlname)
    tuple path(kgenome)        , path(kgenome_i)
    tuple path(dbsnpindel)     , path(dbsnpindel_i)
    tuple path(exac)           , path(exac_i)
    tuple path(evs)            , path(evs_i)
    tuple path(localcontrolwgs), path(localcontrolwgs_i)
    tuple path(localcontrolwes), path(localcontrolwes_i)
    tuple path(gnomadgenomes)  , path(gnomadgenomes_i)
    tuple path(gnomadexomes)   , path(gnomadexomes_i)
    val (chrprefix)

    output:
    tuple val(meta), path('*.ForAnnovar.bed')                         , emit: forannovar
    tuple val(meta), path('*.vcf')                                    , emit: unziped_vcf
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    """
    zcat < $vcf | \\
    annotate_vcf.pl -a - -b $dbsnpindel --columnName='DBSNP' \\
        --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=${params.padding} \\
        --minOverlapFraction=${params.minoverlapfraction} --maxBorderDistanceSum=${params.maxborderdist} \\
        --maxNrOfMatches=${params.maxmatches} | \\
    annotate_vcf.pl -a - -b $kgenome --columnName='1K_GENOMES' \\
        --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=${params.padding} \\
        --minOverlapFraction=${params.minoverlapfraction} --maxBorderDistanceSum=${params.maxborderdist} \\
        --maxNrOfMatches=${params.maxmatches} | \\
    annotate_vcf.pl -a - -b $exac --columnName='ExAC' \\
        --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $evs --columnName='EVS' \\
        --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $gnomadexomes --columnName='GNOMAD_EXOMES' \\
        --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $gnomadgenomes --columnName='GNOMAD_GENOMES' \\
        --bFileType vcf --reportLevel 4 --reportMatchType| \\
    annotate_vcf.pl -a - -b $localcontrolwgs --columnName='LocalControlAF_WGS' \\
        --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType | \\
    annotate_vcf.pl -a - -b $localcontrolwes --columnName='LocalControlAF_WES' \\
        --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType | \\
    tee ${prefix}.vcf | vcf_to_annovar.pl $chrprefix "" > ${prefix}.ForAnnovar.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

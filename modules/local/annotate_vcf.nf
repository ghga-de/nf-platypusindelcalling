//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v1' :'kubran/odcf_platypusindelcalling:v1' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple val(meta2),path(kgenome),path(kgenome_i),path(dbsnpindel),path(dbsnpindel_i),path(exac),path(evs),path(evs_i),path(exac_i),path(localcontrolwgs),path(localcontrolwgs_i),path(localcontrolwes),path(localcontrolwes_i),path(gnomadgenomes),path(gnomadgenomes_i),path(gnomadexomes),path(gnomadexomes_i)
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
    def pipe  = [dbsnpindel ? " | annotate_vcf.pl -a - -b ${dbsnpindel} --columnName='DBSNP' --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=${params.padding} --minOverlapFraction=${params.minoverlapfraction} --maxBorderDistanceSum=${params.maxborderdist} --maxNrOfMatches=${params.maxmatches}" : '',
                kgenome ? " | annotate_vcf.pl -a - -b ${kgenome} --columnName='1K_GENOMES' --reportMatchType --bAdditionalColumn=2 --reportBFeatCoord --padding=${params.padding} --minOverlapFraction=${params.minoverlapfraction} --maxBorderDistanceSum=${params.maxborderdist} --maxNrOfMatches=${params.maxmatches}" : '',
                exac ? " | annotate_vcf.pl -a - -b ${exac} --columnName='ExAC' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                evs ? " | annotate_vcf.pl -a - -b ${evs} --columnName='EVS' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                localcontrolwgs ? " | annotate_vcf.pl -a - -b ${localcontrolwgs} --columnName='LocalControlAF_WGS' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                localcontrolwes ? " | annotate_vcf.pl -a - -b ${localcontrolwes} --columnName='LocalControlAF_WES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadgenomes ? " | annotate_vcf.pl -a - -b ${gnomadgenomes} --columnName='GNOMAD_GENOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadexomes ? " | annotate_vcf.pl -a - -b ${gnomadexomes} --columnName='GNOMAD_EXOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : ''
                ].join(' ').trim()

    """
    zcat < $vcf $pipe |\\
        tee ${prefix}.vcf | vcf_to_annovar.pl $chrprefix "" > ${prefix}.ForAnnovar.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

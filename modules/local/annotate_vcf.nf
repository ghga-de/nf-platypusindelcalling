//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v4.sif' :
    'kubran/odcf_indelcalling:v4' }"

    publishDir params.outdir+'/annotate_vcf' , mode: 'copy'

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)
    tuple path(kgenome), path(kgenome_i)
    tuple path(dbsnpindel), path(dbsnpindel_i)
    tuple path(dbsnpsnv), path(dbsnpsnv_i)
    tuple path(exac), path(exac_i)
    tuple path(evs), path(evs_i)
    tuple path(localcontrolwgs), path(localcontrolwgs_i)
    tuple path(localcontrolwes), path(localcontrolwes_i)
    tuple path(gnomadgenomes), path(gnomadgenomes_i)
    tuple path(gnomadexomes), path(gnomadexomes_i)

    output:
    tuple val(meta), path('*.ForAnnovar.bed')                    , emit: forannovar
    tuple val(meta), path('*.vcf')                               , emit: unziped_vcf
    path  "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def for_annovar = "${meta.id}.ForAnnovar.bed"
    def temp_name   = "${meta.id}.tmp"
    def out_vcf     = "${meta.id}.vcf"

    // QUESTION: BETTER WAY TO RUN ANNOTATE_VCF, either use the parameters check or run as block

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
    tee $temp_name | vcf_to_annovar.pl ${params.chr_prefix} ${params.chr_suffix} > $for_annovar

    mv $temp_name $out_vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    END_VERSIONS
    """
}

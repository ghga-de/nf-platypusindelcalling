//# Gene annotation with annovar
// PROCESS ANNOVAR table_annovar
// working database is annovar_Feb2016

process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"

    publishDir params.outdir+'/annotate_vcf'                           , mode: 'copy'
    
    input:
    tuple val(meta), file(annovar_bed)
    tuple val(meta), file(ch_vcf)
    each file(annovar_table)
    val(chrprefix)

    output:
    tuple val(meta), path('*.temp.vcf.gz'), path('*.temp.vcf.gz.tbi')   , emit: vcf
    tuple val(meta), path('*.log')                                      , emit: log
    path  "versions.yml"                                                , emit: versions
    tuple val(meta), path('*_genomicSuperDups')                         
    tuple val(meta), path('*_cytoBand')                                 
    tuple val(meta), path('*variant_function')                          
    tuple val(meta), path('*exonic_variant_function')                  

    when:
    task.ext.when == null || task.ext.when

    script:
    def av_segdup        = "${meta.id}.${params.buildver}_genomicSuperDups"
    def av_cytoband      = "${meta.id}.${params.buildver}_cytoBand"
    def newcol           = "${meta.id}.newcol.tsv"
    def chr_prefix  = chrprefix == "dummy" ? "" : chrprefix 

    """
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -dbtype ${params.dbtype} $annovar_bed $annovar_table
    ## segdup annotation with annovar
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -regionanno -dbtype segdup --outfile=$meta.id $annovar_bed $annovar_table
    ## cytobad annotation with annovar
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -regionanno -dbtype band --outfile=$meta.id $annovar_bed $annovar_table
    
    processAnnovarOutput.pl "${meta.id}.ForAnnovar.bed.variant_function" "${meta.id}.ForAnnovar.bed.exonic_variant_function" > $newcol

    newCols2vcf.pl --vcfFile=$ch_vcf --newColFile=$newcol \\
        --newColHeader=${params.geneannocols} \\
        --chrPrefix=$chr_prefix \\
        --chrSuffix="" \\
        --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_segdup \\
        --newColHeader=${params.segdupcol} \\
        --chrPrefix=$chr_prefix \\
        --chrSuffix="" \\
        --reportColumns="1" --bChrPosEnd="2,7,8" |
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_cytoband \\
        --newColHeader=${params.cytobandcol} \\
        --chrPrefix=$chr_prefix \\
        --chrSuffix="" \\
        --reportColumns="1" --bChrPosEnd="2,7,8" > "${meta.id}.temp.vcf"

    bgzip -c "${meta.id}.temp.vcf" > "${meta.id}.temp.vcf.gz"
    tabix "${meta.id}.temp.vcf.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

//# Gene annotation with annovar
// PROCESS ANNOVAR table_annovar
// working database is annovar_Feb2016

process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v4.sif' :
    'kubran/odcf_indelcalling:v4' }"

    publishDir params.outdir+'/annotate_vcf' , mode: 'copy'
    
    input:
    tuple val(meta), file(annovar_bed)
    tuple val(meta), file(ch_vcf)
    each file(annovar_table)

    output:
    tuple val(meta), path('*.temp.vcf.gz'), path('*.temp.vcf.gz.tbi')   , emit: vcf
    tuple val(meta), path('*.log')                                      , emit: log
    tuple val(meta), path('*_genomicSuperDups')                         
    tuple val(meta), path('*_cytoBand')                                 
    tuple val(meta), path('*variant_function')                          , emit: variant_function
    tuple val(meta), path('*exonic_variant_function')                   , emit: exonic_variant_function
    path  "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def variant_function = "${meta.id}.ForAnnovar.bed.variant_function"
    def variant_exon     = "${meta.id}.ForAnnovar.bed.exonic_variant_function"
    def tempname         = "${meta.id}.temp.vcf"
    def tempnamegz       = "${meta.id}.temp.vcf.gz"
    def av_segdup        = "${meta.id}.${params.buildver}_genomicSuperDups"
    def av_cytoband      = "${meta.id}.${params.buildver}_cytoBand"
    def newcol           = "${meta.id}.newcol.tsv" 

    """
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -dbtype ${params.dbtype} $annovar_bed $annovar_table
    ## segdup annotation with annovar
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -regionanno -dbtype segdup --outfile=$meta.id $annovar_bed $annovar_table
    ## cytobad annotation with annovar
    perl ${params.annovar_path}/annotate_variation.pl \\
        --buildver=${params.buildver} -regionanno -dbtype band --outfile=$meta.id $annovar_bed $annovar_table
    
    processAnnovarOutput.pl $variant_function $variant_exon > $newcol

    newCols2vcf.pl --vcfFile=$ch_vcf --newColFile=$newcol \\
        --newColHeader=${params.geneannocols} \\
        --chrPrefix=${params.chr_prefix} \\
        --chrSuffix=${params.chr_suffix} \\
        --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_segdup \\
        --newColHeader=${params.segdupcol} \\
        --chrPrefix=${params.chr_prefix} \\
        --chrSuffix=${params.chr_suffix} \\
        --reportColumns="1" --bChrPosEnd="2,7,8" |
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_cytoband \\
        --newColHeader=${params.cytobandcol} \\
        --chrPrefix=${params.chr_prefix} \\
        --chrSuffix=${params.chr_suffix} \\
        --reportColumns="1" --bChrPosEnd="2,7,8" > $tempname

    bgzip -c $tempname > $tempnamegz
    tabix $tempnamegz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    annovar_Feb2016: ${params.buildver}  
    perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
    gzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^.*gzip //; s/Using.*\$//')
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

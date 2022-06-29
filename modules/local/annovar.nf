// PROCESS ANNOVAR table_annovar
// working database is annovar_Feb2016

process ANNOVAR {
    tag "$sample"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '' :
    'kubran/indelcalling_perl:v4' }"

//conda wont work here, database files are embeded inside the docker file
//'noelnamai/annovar:4.18' is also a docker from https://github.com/noelnamai/aws-mutation-calling/blob/master/main.nf

    publishDir params.outdir+'/annovar' , mode: 'copy'

    input:
    tuple val(sample), file(annovar_bed)
    tuple val(sample), file(ch_vcf)
    each file(annovar_table)


    output:
    tuple val(sample), path('*.temp.vcf.gz'), path('*.temp.vcf.gz.tbi')   , emit: vcf
    tuple val(sample), path('*.log')                                      , emit: log
    tuple val(sample), path('*_genomicSuperDups')                         , emit: segdup
    tuple val(sample), path('*_cytoBand')                                 , emit: cytobad
    tuple val(sample), path('*variant_function')                          , emit: variant_function
    tuple val(sample), path('*exonic_variant_function')                   , emit: exonic_variant_function
//    path  "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def variant_function ="${sample}.ForAnnovar.bed.variant_function"
    def variant_exon="${sample}.ForAnnovar.bed.exonic_variant_function"
    def tempname = "${sample}.temp.vcf"
    def tempnamegz = "${sample}.temp.vcf.gz"
    def av_segdup =  "${sample}.${params.buildver}_genomicSuperDups"
    def av_cytoband = "${sample}.${params.buildver}_cytoband"

    """
    perl ${params.annovar_path}/annotate_variation.pl --buildver=${params.buildver} -dbtype ${params.dbtype} $annovar_bed $annovar_table
    perl ${params.annovar_path}/annotate_variation.pl --buildver=${params.buildver} -regionanno -dbtype segdup --outfile=$sample $annovar_bed $annovar_table
    perl ${params.annovar_path}/annotate_variation.pl --buildver=${params.buildver} -regionanno -dbtype band --outfile=$sample $annovar_bed $annovar_table

    processAnnovarOutput.pl $variant_function $variant_exon > newcol.tsv

    newCols2vcf.pl --vcfFile=$ch_vcf --newColFile='newcol.tsv' \\
        --newColHeader=${params.geneannocols} --chrPrefix=${params.chr_prefix} --chrSuffix=${params.chr_suffix} --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" | \\
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_segdup --newColHeader=${params.segdupcol} --chrPrefix=${params.chr_prefix} --chrSuffix=${params.chr_suffix} \\
        --reportColumns="1" --bChrPosEnd="2,7,8"  | \\
    newCols2vcf.pl --vcfFile="-" --newColFile=$av_cytoband --newColHeader=${params.cytobandcol} --chrPrefix=${params.chr_prefix} --chrSuffix=${params.chr_suffix} \\
        --reportColumns="1" --bChrPosEnd="2,7,8" > $tempname

    bgzip -c $tempname > $tempnamegz
    tabix $tempnamegz
    """
}

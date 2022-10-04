//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATION_PIPES {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v4.sif' :
    'kubran/odcf_indelcalling:v4' }"

    publishDir params.outdir+'/annotation_pipe' , mode: 'copy'

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
    tuple val(meta), path('*.deepanno.vcf.gz'), path('*.deepanno.vcf.gz.tbi') , emit: deep_vcf
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def enchangers_seq = enchangers ? "-en ${enchangers}" : ''
    def cpgislands_seq = cpgislands ? "-cp ${cpgislands}" : ''
    def tfbscons_seq = tfbscons ? "-tf ${tfbscons}" : ''
    def mirnas_snornas_seq = mirnas_snornas ? "-ms ${mirnas_snornas}" : '' 
    def encode_dnase_seq = encode_dnase ? "-ed ${encode_dnase}" : ''
    def mirbase_seq = mirbase ? "-mir ${mirbase}" : ''
    def cosmic_seq = cosmic ? "-c ${cosmic}" : '' 
    def mir_targets_seq = mir_targets ? "-mt ${mir_targets}" : '' 
    def cgi_mountains_seq = cgi_mountains ? "-cm ${cgi_mountains}" : ''
    def phastconselem_seq = phastconselem ? "-p ${phastconselem}" : ''
    def encode_tfbs_seq = encode_tfbs ? "-et ${encode_tfbs}" : ''


    """
    create_pipes.sh -i $vcf -id ${meta.id} $enchangers_seq $cpgislands_seq $tfbscons_seq $mirnas_snornas_seq $encode_dnase_seq $mirbase_seq $cosmic_seq $mir_targets_seq $cgi_mountains_seq $phastconselem_seq $encode_tfbs_seq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

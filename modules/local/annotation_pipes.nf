//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATION_PIPES {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v7.sif' :'kubran/odcf_indelcalling:v7' }"

    debug true
    
    input:
    tuple val(meta)           , file(vcf)             , file(vcf_tbi)
    tuple file(enchangers)    , file(enchangers_i)
    tuple file(cpgislands)    , file(cpgislands_i)
    tuple file(tfbscons)      , file(tfbscons_i)
    tuple file(encode_dnase)  , file(encode_dnase_i)
    tuple file(mirnas_snornas), file(mirnas_snornas_i)
    tuple file(cosmic)        , file(cosmic_i)
    tuple file(mirbase)       , file(mirbase_i)
    tuple file(mir_targets)   , file(mir_targets_i)
    tuple file(cgi_mountains) , file(cgi_mountains_i)
    tuple file(phastconselem) , file(phastconselem_i)
    tuple file(encode_tfbs)   , file(encode_tfbs_i)

    output:
    tuple val(meta), path('*.deepanno.vcf.gz'), path('*.deepanno.vcf.gz.tbi') , emit: deep_vcf
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def enchangers_seq     = enchangers ? "-en ${enchangers}" : ''
    def cpgislands_seq     = cpgislands ? "-cp ${cpgislands}" : ''
    def tfbscons_seq       = tfbscons ? "-tf ${tfbscons}" : ''
    def mirnas_snornas_seq = mirnas_snornas ? "-ms ${mirnas_snornas}" : '' 
    def encode_dnase_seq   = encode_dnase ? "-ed ${encode_dnase}" : ''
    def mirbase_seq        = mirbase ? "-mir ${mirbase}" : ''
    def cosmic_seq         = cosmic ? "-c ${cosmic}" : '' 
    def mir_targets_seq    = mir_targets ? "-mt ${mir_targets}" : '' 
    def cgi_mountains_seq  = cgi_mountains ? "-cm ${cgi_mountains}" : ''
    def phastconselem_seq  = phastconselem ? "-p ${phastconselem}" : ''
    def encode_tfbs_seq    = encode_tfbs ? "-et ${encode_tfbs}" : ''


    def pipe               = [enchangers_seq,
                            cpgislands_seq,
                            tfbscons_seq,
                            mirnas_snornas_seq,
                            encode_dnase_seq,
                            mirbase_seq,
                            cosmic_seq,
                            mir_targets_seq,
                            cgi_mountains_seq,
                            phastconselem_seq,
                            encode_tfbs_seq].join(' ').trim() 

    """
    create_pipes.sh -i $vcf -id ${prefix} $pipe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

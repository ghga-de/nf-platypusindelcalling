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
    tuple val(meta), path('*.deepanno.vcf.gz'), path('*.deepanno.vcf.gz.tbi') , emit: vcf
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def pipe               = [enchangers.baseName !='input' ? "-en ${enchangers}" : '',
                            cpgislands.baseName !='input' ? "-cp ${cpgislands}" : '',
                            tfbscons.baseName !='input' ? "-tf ${tfbscons}" : '',
                            mirnas_snornas.baseName !='input' ? "-ms ${mirnas_snornas}" : '',
                            encode_dnase.baseName !='input' ? "-ed ${encode_dnase}" : '',
                            mirbase.baseName !='input' ? "-mir ${mirbase}" : '',
                            cosmic.baseName !='input' ? "-c ${cosmic}" : '',
                            mir_targets.baseName !='input' ? "-mt ${mir_targets}" : '' ,
                            cgi_mountains.baseName !='input' ? "-cm ${cgi_mountains}" : '',
                            phastconselem.baseName !='input' ? "-p ${phastconselem}" : '',
                            encode_tfbs.baseName !='input' ? "-et ${encode_tfbs}" : ''].join(' ').trim() 
    """
    create_pipes.sh -i $vcf -id ${prefix} $pipe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}


process SAMPLE_SWAP {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v1' :'kubran/odcf_platypusindelcalling:v1' }"

    input:
    tuple val(meta) , file(ch_vcf), file(ch_vcf_i),  val(tumorname), val(controlname)
    tuple path(ref) , path(ref_fai)
    each path(chrlength_file)
    each path(genemodel)              
    tuple path(localcontroltindawgs), path(localcontroltindawgs_tbi)
    tuple path(localcontroltindawes), path(localcontroltindawes_tbi)
    tuple path(gnomadgenomes)       , path(gnomadgenomes_tbi)
    tuple path(gnomadexomes)        , path(gnomadexomes_tbi)
    val chrprefix

    output:
    tuple val(meta), path('indel_*.tinda.vcf')                           , emit: vcf  , optional: true
    tuple val(meta), path('indel_*.swap.json')                           , emit: json , optional: true
    path  "snvs_*.GTfiltered_raw.vcf"                                    , optional: true
    path  "snvs_*.GTfiltered_gnomAD.vcf"                                 , optional: true
    path  "snvs_*.GTfiltered_gnomAD.SomaticIn.vcf"                       , optional: true   
    path  "snvs_*.GTfiltered_gnomAD.Germline.Rare.vcf"                   , optional: true
    path  "snvs_*.GTfiltered_gnomAD.Germline.Rare.txt"                   , optional: true
    path  "snvs_*.GTfiltered_gnomAD.Germline.Rare.Rescue.png"            , optional: true
    path  "snvs_*.GTfiltered_gnomAD.Germline.Rare.Rescue.txt"            , optional: true
    path  "indel_*.checkSampleSwap_TiN.log"                              , emit: log
    path  "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    if (meta.iscontrol == '1'){      
        """
        checkSampleSwap_TiN.pl \\
            --pid=$prefix \\
            --raw_file=$ch_vcf \\
            --chr_prefix=$chrprefix \\
            --gnomAD_genome=$gnomadgenomes \\
            --gnomAD_exome=$gnomadexomes \\
            --localControl_WGS=$localcontroltindawgs \\
            --localControl_WES=$localcontroltindawes \\
            --TiNDA_rightBorder=${params.right_border} \\
            --TiNDA_bottomBorder=${params.bottom_border} \\
            --maf_thershold=${params.maf_threshold} \\
            --chrLengthFile=$chrlength_file \\
            --normal_header_col=$controlname \\
            --tumor_header_col=$tumorname \\
            --sequenceType=${params.seqtype} \\
            --gene_model_bed=$genemodel \\
            --reference=$ref \\
            --outfile_tindaVCF=indel_${prefix}.tinda.vcf \\
            --outfile_swapJSON=indel_${prefix}.swap.json \\
            2>&1 | tee indel_${prefix}.checkSampleSwap_TiN.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        touch indel_empty.checkSampleSwap_TiN.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}

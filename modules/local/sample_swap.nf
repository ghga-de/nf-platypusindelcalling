
process SAMPLE_SWAP {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v6.sif' :
    'kubran/odcf_indelcalling:v6' }"

//    publishDir params.outdir+'/sample_swap' , mode: 'copy'

    input:
    tuple val(meta)      , file(ch_vcf), file(ch_vcf_i),  val(tumorname), val(controlname)
    tuple path(ref)                    , path(ref_fai)
    each path(chrlength_file)
    tuple path(genemodel)              , path(genemodel_tbi)
    tuple path(localcontrolplatypuswgs), path(localcontrolplatypuswgs_tbi)
    tuple path(localcontrolplatypuswes), path(localcontrolplatypuswes_tbi)
    tuple path(gnomadgenomes)          , path(gnomadgenomes_tbi)
    tuple path(gnomadexomes)           , path(gnomadexomes_tbi)
    val chrprefix

    output:
    tuple val(meta), path('indel_*.tinda.vcf')                           , optional: true
    tuple val(meta), path('indel_*.swap.json')                           , optional: true
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
    def chr_prefix  = chrprefix == "dummy" ? " " : chrprefix
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    if (meta.iscontrol == '1'){
        
        """
        checkSampleSwap_TiN.pl \\
        --pid=$prefix \\
        --raw_file=$ch_vcf \\
        --chr_prefix=$chr_prefix \\
        --gnomAD_genome=$gnomadgenomes \\
        --gnomAD_exome=$gnomadexomes \\
        --localControl_WGS=$localcontrolplatypuswgs \\
        --localControl_WES=$localcontrolplatypuswes \\
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
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        touch indel_${prefix}.checkSampleSwap_TiN.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}

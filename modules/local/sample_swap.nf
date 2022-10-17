
process SAMPLE_SWAP {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"

    publishDir params.outdir+'/sample_swap' , mode: 'copy'

    input:
    tuple val(meta), file(ch_vcf), file(ch_vcf_i)
    tuple path(ref), path(ref_fai)
    each path(chrlength_file)
    tuple path(genemodel), path(genemodel_tbi)
    tuple path(localcontrolplatypuswgs), path(localcontrolplatypuswgs_tbi)
    tuple path(localcontrolplatypuswes), path(localcontrolplatypuswes_tbi)
    tuple path(gnomadgenomes), path(gnomadgenomes_tbi)
    tuple path(gnomadexomes), path(gnomadexomes_tbi)
    val (chrprefix)

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
    def chr_prefix  = chrprefix == "dummy" ? "" : chrprefix

    if (meta.iscontrol == '1'){
        
        """
        samtools view -H $meta.control_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.controlname.txt
        samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.tumorname.txt

        checkSampleSwap_TiN.pl \\
        --pid=$meta.id \\
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
        --normal_header_col=\$(cat ${meta.id}.controlname.txt) \\
        --tumor_header_col=\$(cat ${meta.id}.tumorname.txt) \\
        --sequenceType=${params.seqtype} \\
        --gene_model_bed=$genemodel \\
        --reference=$ref \\
        --outfile_tindaVCF=indel_${meta.id}.tinda.vcf \\
        --outfile_swapJSON=indel_${meta.id}.swap.json \\
        2>&1 | tee indel_${meta.id}.checkSampleSwap_TiN.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        touch indel_${meta.id}.checkSampleSwap_TiN.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}

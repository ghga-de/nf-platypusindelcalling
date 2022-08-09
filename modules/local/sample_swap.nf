
process SAMPLE_SWAP {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "$baseDir/assets/perlenvironment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'library://kubran/indelcalling/odcf_indelcalling:v0' :
    'kubran/odcf_indelcalling:v0' }"

    publishDir params.outdir+'/tinda'                              , mode: 'copy'

    input:
    tuple val(meta), file(ch_vcf), file(ch_vcf_i)
    tuple path(ref), path(ref_fai)
    each path(chrlength)
    tuple path(genemodel), path(genemodel_tbi)
    tuple path(localcontrolwgs), path(localcontrolwgs_tbi)
    tuple path(localcontrolwes), path(localcontrolwes_tbi)
    tuple path(gnomadgenomes), path(gnomadgenomes_tbi)
    tuple path(gnomadexomes), path(gnomadexomes_tbi)

    output:
    tuple val(meta), path('*.tinda.vcf')                           , emit: vcf
    tuple val(meta), path('*.swap.json')                           , emit: json
    tuple val(meta), path('*.checkSampleSwap_TiN.log')             , emit: log
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta.iscontrol == '1'){
        """
        samtools view -H $meta.control_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.controlname.txt
        samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq > ${meta.id}.tumorname.txt

        checkSampleSwap_TiN.pl --pid=$meta.id --raw_file=$ch_vcf --chr_prefix=${params.chr_prefix} --gnomAD_genome=$gnomadgenomes \\
            --gnomAD_exome=$gnomadexomes --localControl_WGS=$localcontrolwgs --localControl_WES=$localcontrolwes \\
            --TiNDA_rightBorder=${params.right_border} --TiNDA_bottomBorder=${params.bottom_border} --maf_thershold=${params.maf_threshold} \\
            --chrLengthFile=$chrlength --normal_header_col=\$(cat ${meta.id}.controlname.txt) --tumor_header_col=\$(cat ${meta.id}.tumorname.txt) \\
            --sequenceType=${params.seqtype} --gene_model_bed=$genemodel --reference $ref\\
            --outfile_tindaVCF=${meta.id}.tinda.vcf --outfile_swapJSON=${meta.id}.swap.json \\
            2>&1 | tee ${meta.id}.checkSampleSwap_TiN.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
        bedtools: \$(echo \$(bedtools --versions 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
    """
    touch ${meta.id}.swap.json
    touch ${meta.id}.checkSampleSwap_TiN.log
    touch ${meta.id}.tinda.vcf
    touch 'versions.yml'
    """
    }
}

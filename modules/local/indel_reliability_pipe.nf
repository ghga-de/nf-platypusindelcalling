//indel_reliability_pipe

process INDEL_RELIABILITY_PIPE {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"
    
    input:
    tuple val(meta),                 file(ch_vcf),               file(ch_vcf_i)
    tuple file(repeatmasker),        file(repeatmasker_i)
    tuple file(dacblacklist),        file(dacblacklist_i)
    tuple file(dukeexcluded),        file(dukeexcluded_i)
    tuple file(hiseqdepth),          file(hiseqdepth_i)
    tuple file(selfchain),           file(selfchain_i)
    tuple file(mapability),          file(mapability_i)
    tuple file(simpletandemrepeats), file(simpletandemrepeats_i)

    output:
    tuple val(meta), path('*.annotated.vcf.gz'), path('*.annotated.vcf.gz.tbi')   , emit: vcf
    path  "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def pipe  = [repeatmasker.baseName !='input' ? " | annotate_vcf.pl -a - -b ${repeatmasker} --bFileType=bed --columnName='REPEAT_MASKER'" : '',
                dacblacklist.baseName !='input' ? " | annotate_vcf.pl -a - -b ${dacblacklist}  --bFileType=bed --columnName='DAC_BLACKLIST'" : '',
                dukeexcluded.baseName !='input' ? " | annotate_vcf.pl -a - -b ${dukeexcluded} --bFileType=bed --columnName='DUKE_EXCLUDED'" : '',
                hiseqdepth.baseName !='input' ? " | annotate_vcf.pl -a - -b ${hiseqdepth} --bFileType=bed --columnName='HISEQDEPTH'" : '',
                selfchain.baseName !='input' ? " | annotate_vcf.pl -a - -b ${selfchain} --bFileType=bed --columnName='SELFCHAIN' --bAdditionalColumns=4 --maxNrOfMatches=5" : '',
                mapability.baseName !='input' ? " | annotate_vcf.pl -a - -b ${mapability} --bFileType=bed --columnName='MAPABILITY' --breakPointMode --aEndOffset=1" : '',
                simpletandemrepeats.baseName !='input' ? " | annotate_vcf.pl -a - -b ${simpletandemrepeats} --bFileType=bed --columnName='SIMPLE_TANDEMREPEATS' --bAdditionalColumns=4" : ''
                ].join(' ').trim()
    """
    zcat < $ch_vcf $pipe > ${prefix}.annotated.vcf

    bgzip -c ${prefix}.annotated.vcf > ${prefix}.annotated.vcf.gz
    tabix ${prefix}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

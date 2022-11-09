// Takes screenshots of somatic functional variants given window sizes
process VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v7.sif' :'kubran/odcf_indelcalling:v7' }"
    
    errorStrategy 'ignore'

    input:
    tuple val(meta)          ,file(vcf)
    tuple path(ref)          ,path(ref_fai)
    tuple path(repeatmasker) ,path(repeatmasker_tbi)

    output:
    tuple val(meta)          , path('*.indel_somatic_functional_combined.pdf')  , optional: true
    path('versions.yml')                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def control = meta.iscontrol == '1' ? 'true': 'false'

    """
    check_variants_size.sh \\
        -i $vcf \\
        -v ${params.max_var_screenshots} \\
        -c $meta.control_bam \\
        -t $meta.tumor_bam \\
        -r $ref \\
        -w ${params.window_size} \\
        -m $repeatmasker \\
        -s $control \\
        -o ${prefix}.indel_somatic_functional_combined.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        Ghostscript: \$(echo \$(gs -v 2>&1) | sed 's/^.*GPL Ghostscript //; s/ .*\$//')
    END_VERSIONS
    """

}

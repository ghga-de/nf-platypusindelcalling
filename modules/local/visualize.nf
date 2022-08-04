// Takes screenshots of somatic functional variants given window sizes
process VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'python2.7.sif' :
    'kubran/odcf_indelcalling:v0' }"

    publishDir params.outdir+'/screenshots' , mode: 'copy'

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple path(ref), path(ref_fai)
    tuple path(repeatmasker), path(repeatmasker_tbi)

    output:
    tuple val(meta), path('*.pdf')                   , optional: true
    path('versions.yml')                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def nocontrol = meta.iscontrol == '1' ? 'false': 'true'

    """
    check_variants_size.sh \\
    -i $vcfgz \\
    -v ${params.max_var_screenshots} \\
    -c $meta.control_bam \\
    -t $meta.tumor_bam \\
    -r $ref \\
    -w ${params.window_size} \\
    -m $repeatmasker \\
    -s $nocontrol \\
    -o ${meta.id}.indel_somatic_functional_combined.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
    END_VERSIONS
    """

}

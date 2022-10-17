process CHECK_IF_CORRUPTED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"

    input:
    tuple val(meta), file(vcf), file(vcf_tbi)

    output:
    path("*.linesCorrupt")                         , optional: true
    path("*.11")                                   , optional: true
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args?: ''
    def nocontrol = meta.iscontrol == '1' ? 'false': 'true'
    
// If there is no control (iscontrol=1), isNoControlWorkflow is false

    """
    corrupted.sh -i $vcf -c $nocontrol
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools -v 2>&1) | sed 's/^.*bcftools: //; s/ .*\$//')
    END_VERSIONS
    """
    
}

// Applied only for no-control cases
process FILTER_BY_CRIT {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v4.sif' :
    'kubran/odcf_indelcalling:v4' }"

    publishDir params.outdir+'/filtered_vcf' , mode: 'copy'

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)

    output:
     tuple val(meta), path('*Filtered.vcf.gz'),  path('*Filtered.vcf.gz.tbi')   , emit: vcf
     path  "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

// Filter variants only if there is no control, else do noting
    if (meta.iscontrol == '1') {
        """
        mv $vcfgz ${meta.id}_noFiltered.vcf.gz
        tabix ${meta.id}_noFiltered.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        if (params.filter_values != "") 
            {
            """
            vcf_filter_by_crit.py $vcfgz ${meta.id}_postFiltered.vcf ${params.filter_values}
            bgzip -c ${meta.id}_postFiltered.vcf > ${meta.id}_postFiltered.vcf.gz
            tabix ${meta.id}_postFiltered.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/Using.*\$//')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
            END_VERSIONS
            """
             }
        else
        {
            """
            mv $vcfgz ${meta.id}_noFiltered.vcf.gz
            tabix ${meta.id}_noFiltered.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*tabix //; s/Using.*\$//')
            END_VERSIONS
            """ 
        }

    }

}

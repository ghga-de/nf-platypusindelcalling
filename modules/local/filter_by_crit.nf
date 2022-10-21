// Applied only for no-control cases
process FILTER_BY_CRIT {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"

 //    publishDir params.outdir+'/filter_by_crit' , mode: 'copy'

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)

    output:
     tuple val(meta), path('*Filtered.vcf.gz'),  path('*Filtered.vcf.gz.tbi')   , emit: vcf
     path  "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"

// Filter variants only if there is no control, else do noting
    if (meta.iscontrol == '1') {
        """
        mv $vcfgz ${prefix}_noFiltered.vcf.gz
        tabix ${prefix}_noFiltered.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
        END_VERSIONS
        """
    }
    else {
        if (params.filter_values != "") 
            {
            """
            vcf_filter_by_crit.py $vcfgz ${prefix}_postFiltered.vcf ${params.filter_values}
            bgzip -c ${prefix}_postFiltered.vcf > ${prefix}_postFiltered.vcf.gz
            tabix ${prefix}_postFiltered.vcf.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
            END_VERSIONS
            """
             }
        else
        {
            """
            mv $vcfgz ${prefix}_noFiltered.vcf.gz
            tabix ${prefix}_noFiltered.vcf.gz
            
            cat <<-END_VERSIONS > versions.yml
             "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
            END_VERSIONS
            """ 
        }

    }

}
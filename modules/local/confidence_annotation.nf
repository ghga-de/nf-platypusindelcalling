process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v1' :'kubran/odcf_platypusindelcalling:v1' }"
    
    input:
    tuple val(meta), val(tumorname), val(controlname), path(vcfgz), path(vcf_tbi)
    val(ref_type)

    output:
    tuple val(meta), path('*.confidence.vcf.gz') ,path('*.confidence.vcf.gz.tbi') , emit: vcf_ann
    path  "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def samples    = meta.iscontrol == "1" ? "--controlColName=$controlname --tumorColName=$tumorname" : "--nocontrol --tumorColName=$tumorname"
    def ref_spec   = ref_type == "hg38" ? "--refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz" : ""
    
    """
    confidenceAnnotation_Indels.py --infile=$vcfgz \\
        $samples \\
        $ref_spec \\
        --gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} \\
        --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} \\
        --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} \\
        --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} \\
        --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf} \\
        | tee indel_${prefix}.confidence.vcf \\
        | cut -f 1-11 > indel_${prefix}.temp.vcf

    bgzip -c indel_${prefix}.confidence.vcf > indel_${prefix}.confidence.vcf.gz
    tabix indel_${prefix}.confidence.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """
    
}

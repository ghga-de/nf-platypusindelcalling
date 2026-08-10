process PREPARE_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_platypusindelcalling:v1' :'kubran/odcf_platypusindelcalling:v1' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control), path(control_bai)
    tuple val(meta2), path(contig_file)
    tuple path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("*_exclude_contigs.bed") , emit: contigs
    path "versions.yml"                            , emit: versions

    script: 
    def reference_flag = tumor.extension == "cram" ? "-T ${fasta}" : ""
    def outfile = "${meta.id}_exclude_contigs.bed"
    
    if (params.contig_file) {
        """
        touch $outfile
        if [ -n "\$(samtools view ${reference_flag} -H $tumor | grep -P "SN:")" ]; then
            samtools view ${reference_flag} -H $tumor | grep -P "SN:" | sed -e 's/@SQ\\tSN://' | cut -f 1 | grep -v -x -P '(chr)?([1-9]|1[0-9]|2[0-2]|[XYM]|MT)' | grep -v -w -F -f <(awk '{print \$1}' $contig_file) > $outfile || true
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """
    }
    else {
        if (params.runcontigs == 'NONE') {
            """
            touch $outfile
            if [ -n "\$(samtools view ${reference_flag} -H $tumor | grep -P "SN:")" ]; then
                samtools view ${reference_flag} -H $tumor | grep -P "SN:" | sed -e 's/@SQ\\tSN://' | cut -f 1 | grep -v -x -P '(chr)?([1-9]|1[0-9]|2[0-2]|[XYM]|MT)' > $outfile || true
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
            END_VERSIONS
            """
        }
        else if (params.runcontigs == 'ALT_HLA') {
            """
            touch $outfile
            if [ -n "\$(samtools view ${reference_flag} -H $tumor | grep -P "SN:")" ]; then
                samtools view ${reference_flag} -H $tumor | grep -P "SN:" | sed -e 's/@SQ\\tSN://' | cut -f 1 | grep -v -x -P '(chr)?([1-9]|1[0-9]|2[0-2]|[XYM]|MT)' | grep -v -P '_alt|HLA' > $outfile || true
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
            END_VERSIONS
            """
        }
        else {
            """
            touch $outfile

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
            END_VERSIONS
            """
        }
    }
}
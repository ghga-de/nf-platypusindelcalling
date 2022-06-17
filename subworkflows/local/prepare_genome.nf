/*
* Prepare genome for benchmarking
*/

params.genome_options   = [:]

include { SAMTOOLS_FAIDX  } from '../../modules/nf-core/modules/samtools/faidx/main'       addParams( options: params.genome_options )

workflow PREPARE_GENOME {
    take:
    ch_fasta

    main:

    /*
    * Index reference if needed
    */
    SAMTOOLS_FAIDX ( ch_fasta )
    ch_fai = SAMTOOLS_FAIDX.out.fai
    samtools_version = SAMTOOLS_FAIDX.out.versions

    ch_fasta
    .join( ch_fai )
    .set{ ch_fasta_fai }

    emit:
    ch_fasta_fai
    samtools_version
}

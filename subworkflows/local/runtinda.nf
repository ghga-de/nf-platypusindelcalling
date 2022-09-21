//
// RUNTINDA: WILL RUN checkSampleSwap_TIN.sh & checkSampleSwap_TIN.pl
//

params.options = [:]

include { SAMPLE_SWAP                                } from '../../modules/local/sample_swap.nf'     addParams( options: params.options )

workflow RUNTINDA {
    take:
    vcf_ch           // channel: [val(meta), vcfgz, vcfgz_tbi]
    ref              // reference channel [ref.fa, ref.fa.fai]
    chrlength        // channel: [file.txt or file.tab]
    genemodel        // channel: [file.bed.gz, file.bed.gz.tbi]
    localcontrolwgs  // channel: [file.bed.gz, file.bed.gz.tbi]
    localcontrolwes  // channel: [file.bed.gz, file.bed.gz.tbi]
    gnomadgenomes    // channel: [file.bed.gz, file.bed.gz.tbi]
    gnomadexomes     // channel: [file.bed.gz, file.bed.gz.tbi]

    main:


    SAMPLE_SWAP(
    vcf_ch, ref, chrlength, genemodel, localcontrolwgs, localcontrolwes, gnomadgenomes, gnomadexomes
    )
    versions = SAMPLE_SWAP.out.versions
    logs = SAMPLE_SWAP.out.log

    emit:
    versions
    logs
    
}

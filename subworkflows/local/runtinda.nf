//
// RUNTINDA: WILL RUN checkSampleSwap_TIN.sh & checkSampleSwap_TIN.pl
//

params.options = [:]

include { SAMPLE_SWAP  } from '../../modules/local/sample_swap.nf'     addParams( options: params.options )

workflow RUNTINDA {
    take:
    vcf_ch                    // channel: [val(meta), vcfgz, vcfgz_tbi, tumorname, controlname]
    ref                       // reference channel [ref.fa, ref.fa.fai]
    chrlength                 // channel: [file.txt or file.tab]
    genemodel                 // channel: [file.bed.gz, file.bed.gz.tbi]
    localcontrolplatypuswgs   // channel: [file.bed.gz, file.bed.gz.tbi]
    localcontrolplatypuswes   // channel: [file.bed.gz, file.bed.gz.tbi]
    gnomadgenomes             // channel: [file.bed.gz, file.bed.gz.tbi]
    gnomadexomes              // channel: [file.bed.gz, file.bed.gz.tbi]
    chrprefix                 // channel: [prefix]          

    main:


    SAMPLE_SWAP(
    vcf_ch, ref, chrlength, genemodel, localcontrolplatypuswgs, localcontrolplatypuswes, 
    gnomadgenomes, gnomadexomes, chrprefix
    )
    versions = SAMPLE_SWAP.out.versions
    logs = SAMPLE_SWAP.out.log

    emit:
    versions
    logs
    
}

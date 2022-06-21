//
// PLATYPUS_RUNNER: RUN PLATYPUS
//

params.options = [:]

include { PLATYPUS                   } from '../../modules/local/platypus.nf'                                 addParams( options: params.options )
include {  TABIX_BGZIPTABIX          } from '../../modules/local/tabix_bgziptabix.nf'          addParams( options: params.options )

workflow PLATYPUS_RUNNER {
    take:
    sample_ch

    main:

    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }

    PLATYPUS (
    sample_ch, ref
    )
    ch_platypus_log = PLATYPUS.out.platypus_log
    platypus_version = PLATYPUS.out.versions

// TODO: CHECK IF VCF FILES CORRUPTED!! This is in indelcalling.sh  it would be different for with or withoutcontrol cases

    TABIX_BGZIPTABIX (
    PLATYPUS.out.platypus_vcf
    )
    ch_platypus_vcf_to_filter_gz = TABIX_BGZIPTABIX.out.gz_tbi
    bgzip_version = TABIX_BGZIPTABIX.out.versions


    emit:
    ch_platypus_vcf_to_filter_gz
    ch_platypus_log
    platypus_version
    bgzip_version
}

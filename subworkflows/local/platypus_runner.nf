//
// PLATYPUS_RUNNER: RUN PLATYPUS
//

params.options = [:]

include { PLATYPUS          } from '../../modules/local/platypus.nf'               addParams( options: params.options )
include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'     addParams( options: params.options )
include { TABIX_BGZIPTABIX  } from '../../modules/local/tabix_bgziptabix.nf'       addParams( options: params.options )

workflow PLATYPUS_RUNNER {
    take:
    sample_ch

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }

    PLATYPUS (
    sample_ch, ref
    )
    vcf_ch=PLATYPUS.out.vcf
    ch_platypus_log = PLATYPUS.out.platypus_log
    platypus_version = PLATYPUS.out.versions

    TABIX_BGZIPTABIX (
    vcf_ch
    )
    ch_vcf = TABIX_BGZIPTABIX.out.gz_tbi
    bgzip_version = TABIX_BGZIPTABIX.out.versions

    CHECK_IF_CORRUPTED (
    ch_vcf
    )

    emit:
    ch_vcf
    ch_platypus_log
    platypus_version
    bgzip_version
}

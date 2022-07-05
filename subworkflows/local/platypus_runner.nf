//
// PLATYPUS_RUNNER: RUN PLATYPUS
//

params.options = [:]

include { PLATYPUS_WITHCONTROL      } from '../../modules/local/platypus_withcontrol.nf'                addParams( options: params.options )
include { PLATYPUS_NOCONTROL        } from '../../modules/local/platypus_NOcontrol.nf'                  addParams( options: params.options )
//include { CHECK_IF_CORRUPTED as CHECK_WITHCONTROL} from '../../modules/local/check_if_corrupted.nf'   addParams( options: params.options )
//include { CHECK_IF_CORRUPTED as CHECK_NOCONTROL} from '../../modules/local/check_if_corrupted.nf'     addParams( options: params.options )
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX1  } from '../../modules/local/tabix_bgziptabix.nf'       addParams( options: params.options )
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX2  } from '../../modules/local/tabix_bgziptabix.nf'       addParams( options: params.options )

workflow PLATYPUS_RUNNER {
    take:
    sample_ch_withcontrol
    sample_ch_nocontrol

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }

    withcontrol_vcf_ch = Channel.empty()
    nocontrol_vcf_ch = Channel.empty()
    ch_platypus_log = Channel.empty()
    platypus_version = Channel.empty()

///// PIPELINE RUN FOR WITH CONTROL SAMPLES /////
    if (sample_ch_withcontrol) {
        PLATYPUS_WITHCONTROL (
        sample_ch_withcontrol, ref
        )
        withcontrol_vcf_ch=PLATYPUS_WITHCONTROL.out.withcontrol_vcf
        ch_platypus_log = ch_platypus_log.mix(PLATYPUS_WITHCONTROL.out.platypus_log)
        platypus_version = platypus_version.mix(PLATYPUS_WITHCONTROL.out.versions)

        TABIX_BGZIPTABIX1 (
        withcontrol_vcf_ch
        )
        ch_vcf_withcontrol = TABIX_BGZIPTABIX1.out.gz_tbi
        bgzip_version = TABIX_BGZIPTABIX1.out.versions
    }
///// PIPELINE RUN FOR WITH NO CONTROL SAMPLES /////
    if (sample_ch_nocontrol) {
        PLATYPUS_NOCONTROL (
        sample_ch_nocontrol, ref
        )
        nocontrol_vcf_ch=PLATYPUS_NOCONTROL.out.nocontrol_vcf
        ch_platypus_log = ch_platypus_log.mix(PLATYPUS_NOCONTROL.out.platypus_log)
        platypus_version = platypus_version.mix(PLATYPUS_NOCONTROL.out.versions)

        TABIX_BGZIPTABIX2 (
        nocontrol_vcf_ch
        )
        ch_vcf_nocontrol = TABIX_BGZIPTABIX2.out.gz_tbi
        bgzip_version = TABIX_BGZIPTABIX2.out.versions
    }


    emit:
    ch_vcf_withcontrol
    ch_vcf_nocontrol
    ch_platypus_log
    platypus_version
    bgzip_version
}

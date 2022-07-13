//
// RUNTILDA: WILL RUN checkSampleSwap_TIN.sh & checkSampleSwap_TIN.pl
//

params.options = [:]

include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'     addParams( options: params.options )
include { TABIX_BGZIPTABIX  } from '../../modules/local/tabix_bgziptabix.nf'       addParams( options: params.options )


workflow RUNTILDA {
    take:
    sample_ch // channel: [val(meta), [tumor bams], [control bams]]

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }

    emit:

}

//
// FILTER VCF: ...
//

params.options = [:]

include { FILTER_BY_CRIT  } from '../../modules/local/filter_by_crit.nf'     addParams( options: params.options )

workflow FILTER_VCF {
    take:
    vcf_ch
    sample_ch

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }

    // RUN vcf_filter_bycrit.pl ONLY IF THERE IS NO CONTROL!
    FILTER_BY_CRIT(
    vcf_ch
    )

    emit:
    vcf_ch
}

//
// FILTER VCF: ...
//

params.options = [:]

include { FILTER_BY_CRIT    } from '../../modules/local/filter_by_crit.nf'       addParams( options: params.options )
include { INDEL_EXTRACTION  } from '../../modules/local/indel_extraction.nf'     addParams( options: params.options )
include { VISUALIZE         } from '../../modules/local/visualize.nf'            addParams( options: params.options )

workflow FILTER_VCF {
    take:
    vcf_ch

    main:
    if (params.reference) { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { ref = Channel.empty() }
    if (params.repeat_masker) { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.empty() }

    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
    vcf_ch
    )

    INDEL_EXTRACTION(
    FILTER_BY_CRIT.out.vcf
    )
    somatic_vcf_ch=INDEL_EXTRACTION.out.vcf

    // CHECK IF THERE IS FUNCTIONAL SOMATIC VARIANTS, IF THERE IS DO FOLLOWING

    VISUALIZE (
    somatic_vcf_ch, ref, repeatmasker
    )

    emit:
    somatic_vcf_ch
}

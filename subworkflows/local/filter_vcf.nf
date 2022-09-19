//
// FILTER VCF: Filtering only applied if the tumor file has no control!
//

params.options = [:]

include { FILTER_BY_CRIT       } from '../../modules/local/filter_by_crit.nf'       addParams( options: params.options )
include { INDEL_EXTRACTION     } from '../../modules/local/indel_extraction.nf'     addParams( options: params.options )
include { VISUALIZE            } from '../../modules/local/visualize.nf'            addParams( options: params.options )
include { INDEL_JSON           } from '../../modules/local/indel_json_v1.nf'        addParams( options: params.options )

workflow FILTER_VCF {
    take:
    vcf_ch        // channel: [val(meta), vcf]
    ref           // reference channel [ref.fa, ref.fa.fai]
    repeatmasker  // channel: [file.bed.gz, file.bed.gz.tbi]


    main:

    versions=Channel.empty()

    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
    vcf_ch
    )
    versions = versions.mix(FILTER_BY_CRIT.out.versions)

    INDEL_EXTRACTION(
    FILTER_BY_CRIT.out.vcf
    )
    versions = versions.mix(INDEL_EXTRACTION.out.versions)

    // CHECK IF THERE IS FUNCTIONAL SOMATIC VARIANTS, IF THERE IS DO FOLLOWING

    VISUALIZE (
    INDEL_EXTRACTION.out.somatic_functional, ref, repeatmasker
    )
    versions = versions.mix(VISUALIZE.out.versions)

    INDEL_JSON(
    INDEL_EXTRACTION.out.somatic_indel
    )
    versions = versions.mix(INDEL_JSON.out.versions)

    emit:
    versions
}

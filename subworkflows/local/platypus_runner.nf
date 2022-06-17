//
// PLATYPUS_RUNNER: RUN PLATYPUS
//

params.options = [:]

include { PLATYPUS             } from '../../modules/local/platypus.nf'             addParams( options: params.options )
//include { BGZIP_TABIX          } from '../../modules/local/bgzip_tabix'          addParams( options: params.options )
//include { BCFTOOLS_FILTER      } from '../../modules/local/bcftools_filter'      addParams( options: params.options )

workflow PLATYPUS_RUNNER {
    take:
    sample_ch

    main:

    PLATYPUS (
        sample_ch, params.reference, params.reference+'.fai'
)
    ch_platypus_vcf_to_filter = PLATYPUS.out.platypus_vcf
    platypus_version = PLATYPUS.out.versions

    emit:
    ch_platypus_vcf_to_filter
    platypus_version

}

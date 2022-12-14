//
// INDEL_CALLING: RUN PLATYPUS
//

params.options = [:]
 
include { PLATYPUS          } from '../../modules/local/platypus.nf'                      addParams( options: params.options )
include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'            addParams( options: params.options )
include { BCFTOOLS_SORT     } from '../../modules/nf-core/modules/bcftools/sort/main'     addParams( options: params.options )


workflow INDEL_CALLING {
    take:
    sample_ch // channel: [val(meta), tumor,tumor_bai, control, control_bai]
    ref       // reference channel [ref.fa, ref.fa.fai]
    
    main:

    versions=Channel.empty()
    // RUN platypus : calls variants
    PLATYPUS (
    sample_ch, ref
    )
    vcf_ch=PLATYPUS.out.vcf
    ch_platypus_log = PLATYPUS.out.platypus_log
    versions = versions.mix(PLATYPUS.out.versions)

    BCFTOOLS_SORT(
        vcf_ch
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)
   
    vcf_ch = BCFTOOLS_SORT.out.vcf 
    //check if the VCF has the correct amount of columns. 
    CHECK_IF_CORRUPTED (
        vcf_ch
    )

    emit:
    vcf_ch
    ch_platypus_log
    versions
}

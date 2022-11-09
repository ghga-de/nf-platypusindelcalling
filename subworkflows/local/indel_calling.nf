//
// INDEL_CALLING: RUN PLATYPUS
//

params.options = [:]
 
include { PLATYPUS          } from '../../modules/local/platypus.nf'               addParams( options: params.options )
include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'     addParams( options: params.options )

workflow INDEL_CALLING {
    take:
    sample_ch // channel: [val(meta), tumor,tumor_bai, control, control_bai]
    ref       // reference channel [ref.fa, ref.fa.fai]
    
    main:

    // RUN platypus : calls variants
    PLATYPUS (
    sample_ch, ref
    )
    vcf_ch=PLATYPUS.out.vcf
    ch_platypus_log = PLATYPUS.out.platypus_log
    platypus_version = PLATYPUS.out.versions

    // check if the VCF has the correct amount of columns. 
    CHECK_IF_CORRUPTED (
    vcf_ch
    )

    emit:
    vcf_ch
    ch_platypus_log
    platypus_version
}

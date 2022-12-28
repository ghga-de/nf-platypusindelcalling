//
// INDEL_CALLING: RUN PLATYPUS
//

params.options = [:]
 
include { PLATYPUS          } from '../../modules/local/platypus.nf'                      addParams( options: params.options )
include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'            addParams( options: params.options )
include { BCFTOOLS_STATS    } from '../../modules/nf-core/modules/bcftools/stats/main'    addParams( options: params.options )

workflow INDEL_CALLING {
    take:
    sample_ch // channel: [val(meta), tumor, tumor_bai, control, control_bai]
    ref       // reference channel [ref.fa, ref.fa.fai]
    
    main:

    versions=Channel.empty()

    // RUN platypus : calls variants
    PLATYPUS (
    sample_ch, ref
    )
    vcf_ch=PLATYPUS.out.vcf
    log_ch = PLATYPUS.out.log
    versions = versions.mix(PLATYPUS.out.versions)

    // Check if VCF has more than 0 variants
    BCFTOOLS_STATS(
        vcf_ch, [], [], []
    )
    versions = versions.mix(BCFTOOLS_STATS.out.versions)
    
    BCFTOOLS_STATS.out.stats
                .join(vcf_ch)
                .filter{meta, stats, vcf -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }  
                .set{ch_vcf}

    
    //check if the VCF has the correct amount of columns. 
    CHECK_IF_CORRUPTED (
        ch_vcf
    )
    vcf_ch=CHECK_IF_CORRUPTED.out.vcf
    versions.mix(CHECK_IF_CORRUPTED.out.versions)

    emit:
    vcf_ch
    log_ch
    versions
}

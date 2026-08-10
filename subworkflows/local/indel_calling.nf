//
// INDEL_CALLING: RUN PLATYPUS
//

params.options = [:]
 
include { PLATYPUS           } from '../../modules/local/platypus.nf'
include { CHECK_IF_CORRUPTED } from '../../modules/local/check_if_corrupted.nf'
include { BCFTOOLS_STATS     } from '../../modules/nf-core/modules/bcftools/stats/main'
include { PREPARE_CONTIGS    } from '../../modules/local/prepare_contigs.nf'
include { FILTER_CONTIGS     } from '../../modules/local/filter_contigs.nf'


workflow INDEL_CALLING {
    take:
    sample_ch // channel: [val(meta), tumor, tumor_bai, control, control_bai]
    ref       // reference channel [ref.fa, ref.fa.fai]
    contigs  // channel: [val(meta), bed]
    
    main:

    versions=Channel.empty()
    //
    // MODULE: PLATYPUS
    //
    // RUN platypus : calls variants
    PLATYPUS (
        sample_ch, 
        ref
    )
    vcf_ch = PLATYPUS.out.vcf
    log_ch = PLATYPUS.out.log
    versions = versions.mix(PLATYPUS.out.versions)

    CHECK_IF_CORRUPTED (
        vcf_ch
    )
    vcf_ch=CHECK_IF_CORRUPTED.out.vcf
    versions.mix(CHECK_IF_CORRUPTED.out.versions)

    /// filter non-standard contigs optionally
    if (params.runcontigs != "ALL") {
        //
        // MODULE: Prepare contigs file if not provided
        //
        PREPARE_CONTIGS(
            sample_ch,
            contigs,
            ref
            )
        versions = versions.mix(PREPARE_CONTIGS.out.versions)

        FILTER_CONTIGS(
            vcf_ch.join(PREPARE_CONTIGS.out.contigs)
        )
        vcf_ch=FILTER_CONTIGS.out.filtered_vcf
        versions = versions.mix(FILTER_CONTIGS.out.versions)
        
    }

    //
    // MODULE: BCFTOOLS STATS
    //
    // Check if VCF has more than 0 variants
    BCFTOOLS_STATS(
        vcf_ch, 
        [], 
        [], 
        []
    )
    versions = versions.mix(BCFTOOLS_STATS.out.versions)
    
    emit:
    vcf_ch
    log_ch
    versions
}

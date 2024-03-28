//
// OUTPUT_STANDARD_VCF: OUTPUT_STANDARD_VCF
//

params.options = [:]

include { CONVERT_TO_VCF as CONVERT_TO_VCF_1 } from '../../modules/local/convert_to_vcf.nf'   addParams( options: params.options )   
include { CONVERT_TO_VCF as CONVERT_TO_VCF_2 } from '../../modules/local/convert_to_vcf.nf'   addParams( options: params.options )   
include { BCFTOOLS_SORT         } from '../../modules/nf-core/modules/bcftools/sort/main'     addParams( options: params.options )          
include { TABIX_TABIX           } from '../../modules/nf-core/modules/tabix/tabix/main'       addParams( options: params.options )             
include { BCFTOOLS_REHEADER     } from '../../modules/nf-core/modules/bcftools/reheader/main' addParams( options: params.options )             

workflow OUTPUT_STANDARD_VCF {
    take:
    vcf_ch // channel: [val(meta), vcf]
    config
    ref    // reference channel [ref.fa, ref.fa.fai]

    
    main:

    versions=Channel.empty()
    //
    // MODULE: CONVERT_TO_VCF
    //
    CONVERT_TO_VCF_1(
        vcf_ch.map{it -> tuple( it[0], it[1], [])},
        config
    )
    versions = versions.mix(CONVERT_TO_VCF_1.out.versions)

    CONVERT_TO_VCF_2(
        CONVERT_TO_VCF_1.out.std_vcf.map{it -> tuple( it[0], it[1], [])},
        config
    )

    ref.map { it -> tuple([id: it[0].baseName], it[1]) }
            .set{fai}

    BCFTOOLS_REHEADER(
        CONVERT_TO_VCF_2.out.std_vcf.map{it -> tuple( it[0], it[1], [], [])},
        fai
    )

    TABIX_TABIX(
        BCFTOOLS_REHEADER.out.vcf
    )
    versions = versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_SORT(
        BCFTOOLS_REHEADER.out.vcf.join(TABIX_TABIX.out.tbi)
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    versions
}

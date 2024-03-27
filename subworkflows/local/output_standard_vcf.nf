//
// OUTPUT_STANDARD_VCF: OUTPUT_STANDARD_VCF
//

params.options = [:]
 
include { CONVERT_TO_VCF        } from '../../modules/local/convert_to_vcf.nf'                addParams( options: params.options )   
include { BCFTOOLS_SORT         } from '../../modules/nf-core/modules/bcftools/sort/main'     addParams( options: params.options )          
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
    CONVERT_TO_VCF(
        vcf_ch,
        config
    )
    versions = versions.mix(CONVERT_TO_VCF.out.versions)

    ref.map { it -> tuple([id: it[0].baseName], it[1]) }
            .set{fai}

    BCFTOOLS_REHEADER(
        CONVERT_TO_VCF.out.std_vcf.map{it -> tuple( it[0], it[1], [], [])},
        fai

    )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions)

    BCFTOOLS_SORT(
        BCFTOOLS_REHEADER.out.vcf
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    versions
}

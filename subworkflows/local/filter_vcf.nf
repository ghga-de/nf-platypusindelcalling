//
// FILTER VCF: Filtering only applied if the tumor file has no control! Extraction and visualizations are performed for all.
//

params.options = [:]

include { FILTER_BY_CRIT       } from '../../modules/local/filter_by_crit.nf'       addParams( options: params.options )
include { INDEL_EXTRACTION     } from '../../modules/local/indel_extraction.nf'     addParams( options: params.options )
include { VISUALIZE            } from '../../modules/local/visualize.nf'            addParams( options: params.options )
include { INDEL_JSON           } from '../../modules/local/indel_json.nf'           addParams( options: params.options )

workflow FILTER_VCF {
    take:
    vcf_ch        // channel: [val(meta), vcf]
    ref           // reference channel [ref.fa, ref.fa.fai]
    repeatmasker  // channel: [file.bed.gz, file.bed.gz.tbi]


    main:

    versions=Channel.empty()

    //
    // MODULE: FILTER_BY_CRIT
    //
    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
        vcf_ch
    )
    versions = versions.mix(FILTER_BY_CRIT.out.versions)

    //
    // MODULE: INDEL_EXTRACTION
    //
    // RUN indel_extractor_v1.pl : extract somatic, funtional and germline variants
    INDEL_EXTRACTION(
        FILTER_BY_CRIT.out.vcf
    )
    versions = versions.mix(INDEL_EXTRACTION.out.versions)

    // filter out the lists if no variant exists for visualization
    INDEL_EXTRACTION.out.somatic_functional
        .filter{meta, somatic_functional -> WorkflowCommons.getNumLinesInFile(somatic_functional) > 1}
        .set{functional_vars}

    INDEL_EXTRACTION.out.somatic_indel
        .filter{meta, somatic_indel -> WorkflowCommons.getNumLinesInFile(somatic_indel) > 1}
        .set{somatic_indel_ch}

    INDEL_EXTRACTION.out.somatic_ncrna
        .filter{meta, somatic_ncrna -> WorkflowCommons.getNumLinesInFile(somatic_ncrna) > 1}
        .set{somatic_ncrna_ch}

    INDEL_EXTRACTION.out.germline_functional
        .filter{meta, germline_functional -> WorkflowCommons.getNumLinesInFile(germline_functional) > 1}
        .set{germline_funct_ch}

    FILTER_BY_CRIT.out.vcf.map{ it -> tuple( it[0], it[1] )}
                .mix(functional_vars)
                .mix(somatic_indel_ch)
                .mix(somatic_ncrna_ch)
                .mix(germline_funct_ch)
                .set{convert_snvs}               
    //
    // MODULE: VISUALIZE
    //
    // RUN: visualize.py : First,checks if there is functional somatic variants to visualize
    // and then checks the size if smaller than the limit creates a pdf wirk screenshots.
    VISUALIZE (
        functional_vars, 
        ref, 
        repeatmasker
    )
    versions = versions.mix(VISUALIZE.out.versions)

    //
    // MODULE: INDEL_JSON
    //
    // RUN: indel_json_v1.0.pl : Prints indel stats
    INDEL_JSON(
        somatic_indel_ch
    )
    versions = versions.mix(INDEL_JSON.out.versions)

    emit:
    convert_snvs
    versions
}
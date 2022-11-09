//
// PLATYPUS INDEL ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF           } from '../../modules/local/annotate_vcf.nf'            addParams( options: params.options )
include { ANNOVAR                } from '../../modules/local/annovar.nf'                 addParams( options: params.options )
include { INDEL_RELIABILITY_PIPE } from '../../modules/local/indel_reliability_pipe.nf'  addParams( options: params.options )
include { CONFIDENCE_ANNOTATION  } from '../../modules/local/confidence_annotation.nf'   addParams( options: params.options )
include { ANNOTATION_PIPES       } from '../../modules/local/annotation_pipes.nf'        addParams( options: params.options )


workflow INDEL_ANNOTATION {
    take:
    vcf_ch               // channel: [val(meta), , vcf.gz, vcf.gz.tbi ,val(tumorname), val(controlname) ]
    kgenome              // channel: [file.vcf.gz, file.vcf.gz.tbi]
    dbsnpindel           // channel: [file.vcf.gz, file.vcf.gz.tbi]
    exac                 // channel: [file.vcf.gz, file.vcf.gz.tbi]
    evs                  // channel: [file.vcf.gz, file.vcf.gz.tbi]
    localcontrolwgs      // channel: [file.vcf.gz, file.vcf.gz.tbi]
    localcontrolwes      // channel: [file.vcf.gz, file.vcf.gz.tbi]
    gnomadgenomes        // channel: [file.vcf.gz, file.vcf.gz.tbi]
    gnomadexomes         // channel: [file.vcf.gz, file.vcf.gz.tbi]
    annodb               // channel: [table_annovar_dir]
    repeatmasker         // channel: [file.bed.gz, file.bed.gz.tbi]
    dacblacklist         // channel: [file.bed.gz, file.bed.gz.tbi]
    dukeexcluded         // channel: [file.bed.gz, file.bed.gz.tbi]
    hiseqdepth           // channel: [file.bed.gz, file.bed.gz.tbi]
    selfchain            // channel: [file.bed.gz, file.bed.gz.tbi]
    mapability           // channel: [file.bed.gz, file.bed.gz.tbi]
    simpletandemrepeats  // channel: [file.bed.gz, file.bed.gz.tbi]
    enchangers           // channel: [file.bed.gz, file.bed.gz.tbi]
    cpgislands           // channel: [file.bed.gz, file.bed.gz.tbi]
    tfbscons             // channel: [file.bed.gz, file.bed.gz.tbi]
    encode_dnase         // channel: [file.bed.gz, file.bed.gz.tbi]
    mirnas_snornas       // channel: [file.bed.gz, file.bed.gz.tbi]
    cosmic               // channel: [file.bed.gz, file.bed.gz.tbi]
    mirbase              // channel: [file.bed.gz, file.bed.gz.tbi]
    mir_targets          // channel: [file.bed.gz, file.bed.gz.tbi]
    cgi_mountains        // channel: [file.bed.gz, file.bed.gz.tbi]
    phastconselem        // channel: [file.bed.gz, file.bed.gz.tbi]
    encode_tfbs          // channel: [file.bed.gz, file.bed.gz.tbi]
    chr_prefix           // val channel: [prefix]

    main:

    versions=Channel.empty()
    logs=Channel.empty() 

    // RUN annotate_vcf.pl: Uses various databases (all mandatory) to annotate variants
    ANNOTATE_VCF (
    vcf_ch, kgenome, dbsnpindel, exac, evs, localcontrolwgs,
    localcontrolwes, gnomadgenomes, gnomadexomes, chr_prefix
    )
    versions  = versions.mix(ANNOTATE_VCF.out.versions)

    // RUN annovar, processAnnovarOutput.pl and newCols2vcf.pl: annovar annotates and classifies the variants, 
    // perl scripts re-creates vcfs. 
    ANNOVAR(
    ANNOTATE_VCF.out.forannovar, ANNOTATE_VCF.out.unziped_vcf, annodb, chr_prefix
    )
    logs     = logs.mix(ANNOVAR.out.log)
    versions = versions.mix(ANNOVAR.out.versions)

    // RUN annotate_vcf.pl : BED files are used to annotate variants
    INDEL_RELIABILITY_PIPE(
    ANNOVAR.out.vcf, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
    )
    versions = versions.mix(INDEL_RELIABILITY_PIPE.out.versions)

    // RUN: confidenceAnnotation_Indels.py : Confidence annotation will be added to the variants
    CONFIDENCE_ANNOTATION(
    INDEL_RELIABILITY_PIPE.out.vcf, vcf_ch
    )
    ann_vcf_ch  = CONFIDENCE_ANNOTATION.out.vcf_ann
    versions    = versions.mix(CONFIDENCE_ANNOTATION.out.versions)

    // RUN annotate_vcf.pl : Uses optional databases to annotate variants, only given databases will be used. 
    if (params.runIndelDeepAnnotation)
    {
        ANNOTATION_PIPES (
        ann_vcf_ch, enchangers, cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets,
        cgi_mountains, phastconselem, encode_tfbs
        )
        ann_vcf_ch  = ANNOTATION_PIPES.out.vcf 
        versions    = versions.mix(ANNOTATION_PIPES.out.versions)

    }
emit:
logs
ann_vcf_ch
versions
}
